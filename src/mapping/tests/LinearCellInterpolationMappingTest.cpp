#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "logging/LogMacros.hpp"
#include "mapping/LinearCellInterpolationMapping.hpp"
#include "mapping/Mapping.hpp"
#include "math/constants.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(LinearCellInterpolationMapping)

BOOST_AUTO_TEST_CASE(Consistent)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;
  using testing::equals;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createData("InDataScalar", 1, 0_dataID);
  int     inDataScalarID = inDataScalar->getID();

  Vertex &inVertexA = inMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &inVertexB = inMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  Vertex &inVertexC = inMesh->createVertex(Eigen::Vector2d(0.0, 1.0));

  inMesh->allocateDataValues();

  Edge &inEdge0 = inMesh->createEdge(inVertexA, inVertexB);
  Edge &inEdge1 = inMesh->createEdge(inVertexB, inVertexC);
  Edge &inEdge2 = inMesh->createEdge(inVertexC, inVertexA);

  inMesh->createTriangle(inEdge0, inEdge1, inEdge2);
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  inValuesScalar << 1.0, 2.0, 3.0;

  BOOST_CHECK(!inMesh->edges().empty());
  /*
  Should pass: triangles are discarded when the mapping doesn't require the in Participant::setMeshTriangle,
  but Mesh::createTriangle doesn't check it. This is not an integration test!
  */

  BOOST_CHECK(!inMesh->triangles().empty());
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1, 2_dataID);
  int     outDataScalarID = outDataScalar->getID();

  // All vertices to test
  // Center of triangle = average
  outMesh->createVertex(Eigen::Vector2d::Constant(1.0 / 3.0));
  // Exact mapping if grid is matching
  outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  outMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  outMesh->createVertex(Eigen::Vector2d(0.0, 1.0));
  // Fallback on NP when slightly outside
  // AB: exact middle (slightly outside and on side of A), BC: 2/3 on side B. CA: check fall-back on edge if out of domain
  outMesh->createVertex(Eigen::Vector2d(0.49, -0.01));
  outMesh->createVertex(Eigen::Vector2d(2.5 / 3, 1. / 3));
  outMesh->createVertex(Eigen::Vector2d(-10.0, 0.25));
  // Check fall back on nearest neighbor
  outMesh->createVertex(Eigen::Vector2d(-0.1, -0.1)); // Fall back to NN expected
  outMesh->createVertex(Eigen::Vector2d(2.5, -1.0));
  outMesh->createVertex(Eigen::Vector2d(2.5, 10.0));

  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::LinearCellInterpolationMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check expected
  Eigen::VectorXd expected(outMesh->nVertices());
  expected << 2.0, 1.0, 2.0, 3.0, 1.49, 2.25, 1.5, 1.0, 2.0, 3.0;
  BOOST_CHECK(equals(expected, outValuesScalar));
}

BOOST_AUTO_TEST_CASE(Conservative)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;
  using testing::equals;

  /* Send forces on some points, map to a triangle and check correct distribution.
     We only apply forces inside the triangle or at its boundary.
     Behavior when projecting is not ideal, see #1304
   */

  // Map to a triangle ABC which represents the lower-left half of a unit square
  const double forceOnMid          = 1.0;
  const double forceOnMidAB        = 2.0;
  const double forceOnA            = 10.0;
  const double forceOnC            = 7.0;
  const double unbalancedforceOnBC = 3.0; // 75% on B, 25% on C
  const double netForce            = forceOnMid + forceOnMidAB + forceOnA + forceOnC + unbalancedforceOnBC;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createData("InDataScalar", 1, 0_dataID);
  int     inDataScalarID = inDataScalar->getID();

  inMesh->createVertex(Eigen::Vector2d(1. / 3, 1. / 3)); // Mid
  inMesh->createVertex(Eigen::Vector2d(0.5, 0.0));       // Mid AB
  inMesh->createVertex(Eigen::Vector2d(0.0, 0.0));       // A
  inMesh->createVertex(Eigen::Vector2d(0.0, 1.0));       // C
  inMesh->createVertex(Eigen::Vector2d(0.75, 0.25));     // Along BC

  inMesh->allocateDataValues();
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  inValuesScalar << forceOnMid, forceOnMidAB, forceOnA, forceOnC, unbalancedforceOnBC;

  // Create output mesh (the triangle ABC)
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1, 2_dataID);
  int     outDataScalarID = outDataScalar->getID();

  Vertex &outVertexA = outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &outVertexB = outMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  Vertex &outVertexC = outMesh->createVertex(Eigen::Vector2d(0.0, 1.0));

  Edge &outEdge0 = outMesh->createEdge(outVertexA, outVertexB);
  Edge &outEdge1 = outMesh->createEdge(outVertexB, outVertexC);
  Edge &outEdge2 = outMesh->createEdge(outVertexC, outVertexA);

  outMesh->createTriangle(outEdge0, outEdge1, outEdge2);
  BOOST_CHECK(!outMesh->edges().empty());
  BOOST_CHECK(!outMesh->triangles().empty());

  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::LinearCellInterpolationMapping mapping(mapping::Mapping::CONSERVATIVE, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check expected
  Eigen::VectorXd expected(outMesh->nVertices());
  const double    expectedA = forceOnMid / 3 + forceOnMidAB / 2 + forceOnA;
  const double    expectedB = forceOnMid / 3 + forceOnMidAB / 2 + unbalancedforceOnBC * 0.75;
  const double    expectedC = forceOnMid / 3 + forceOnC + unbalancedforceOnBC * 0.25;
  expected << expectedA, expectedB, expectedC;
  // Check expected value, and conservation
  BOOST_CHECK(equals(expected, outValuesScalar));
  BOOST_CHECK(equals(netForce, outValuesScalar.sum()));
}

BOOST_AUTO_TEST_CASE(ConsistentOneTetra3D)
{
  PRECICE_TEST(1_rank);
  int dimensions = 3;
  using testing::equals;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createData("InDataScalar", 1, 0_dataID);
  int     inDataScalarID = inDataScalar->getID();

  Vertex &inVertexA = inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &inVertexB = inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  Vertex &inVertexC = inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  Vertex &inVertexD = inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));

  inMesh->allocateDataValues();

  // Create a tetra
  inMesh->createTetrahedron(inVertexA, inVertexB, inVertexC, inVertexD);
  // Create triangle in the plane z = 0
  inMesh->createTriangle(inVertexA, inVertexB, inVertexC);
  // Add edge BD to check fall-back on edge
  inMesh->createEdge(inVertexB, inVertexD);

  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  inValuesScalar << 1.0, 2.0, 3.0, 4.0; //1 + x + 2y + 3z

  BOOST_CHECK(!inMesh->tetrahedra().empty());
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1, 2_dataID);
  int     outDataScalarID = outDataScalar->getID();

  // Center and the 4 vertices (expected: average then 4 values)
  outMesh->createVertex(Eigen::Vector3d::Constant(0.25));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  // Point below the triangle ABC -> fallback to triangle. Expected projection on triangle.
  outMesh->createVertex(Eigen::Vector3d(1.0 / 3, 1.0 / 3, -0.1));
  // Point close to the triangle ACD (which isn't set!).
  // Wanted behavior: NN => 3.0. Incorrect behavior: 3.6 in case of fall-back to edge. See issue #1304
  outMesh->createVertex(Eigen::Vector3d(-0.1, 0.8, 0.5));
  // Point inside the triangle BCD (not set) -> Check it is inside the tetra. Expected: 0.2*2+0.3*3+0.5*4 = 3.3
  outMesh->createVertex(Eigen::Vector3d(0.2, 0.3, 0.5));
  // Point on the the edge BD for fall-back. Expected: 0.4*3 + 0.6*4 = 3.6
  outMesh->createVertex(Eigen::Vector3d(0, 0.4, 0.6));

  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::LinearCellInterpolationMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check expected
  Eigen::VectorXd expected(outMesh->nVertices());
  expected << 2.5, 1.0, 2.0, 3.0, 4.0, 2.0, 3.0, 3.3, 3.6;
  BOOST_CHECK(equals(expected, outValuesScalar));
}

BOOST_AUTO_TEST_CASE(ConservativeOneTetra3D)
{
  PRECICE_TEST(1_rank);
  int dimensions = 3;
  using testing::equals;

  /* Send forces on some points, map to a tetra and check correct distribution.
     We only apply forces inside the tetra or at its boundary.
     Behavior when projecting is not ideal, see #1304
   */

  const double forceOnMid              = 1.0;
  const double forceOnMidAB            = 2.0;
  const double forceOnA                = 10.0;
  const double forceOnC                = 7.0;
  const double unbalancedforceOnBC     = 3.0;  // 75% on B, 25% on C
  const double unbalancedInternalForce = 11.0; // With barycentric coordinates (0.4, 0.1, 0.2, 0.3)
  const double netForce                = forceOnMid + forceOnMidAB + forceOnA + forceOnC + unbalancedforceOnBC + unbalancedInternalForce;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createData("InDataScalar", 1, 0_dataID);
  int     inDataScalarID = inDataScalar->getID();

  inMesh->createVertex(Eigen::Vector3d(0.25, 0.25, 0.25)); // Mid
  inMesh->createVertex(Eigen::Vector3d(0.5, 0.0, 0.0));    // Mid AB
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));    // A
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));    // C
  inMesh->createVertex(Eigen::Vector3d(0.75, 0.25, 0.0));  // Along BC
  inMesh->createVertex(Eigen::Vector3d(0.1, 0.2, 0.3));

  inMesh->allocateDataValues();
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  inValuesScalar << forceOnMid, forceOnMidAB, forceOnA, forceOnC, unbalancedforceOnBC, unbalancedInternalForce;

  // Create output mesh (the tetra ABCD)
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1, 2_dataID);
  int     outDataScalarID = outDataScalar->getID();

  Vertex &outVertexA = outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &outVertexB = outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  Vertex &outVertexC = outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  Vertex &outVertexD = outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));

  outMesh->createTetrahedron(outVertexA, outVertexB, outVertexC, outVertexD);
  BOOST_CHECK(!outMesh->tetrahedra().empty());

  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::LinearCellInterpolationMapping mapping(mapping::Mapping::CONSERVATIVE, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check expected
  Eigen::VectorXd expected(outMesh->nVertices());
  const double    expectedA = forceOnMid / 4 + forceOnMidAB / 2 + forceOnA + 0.4 * unbalancedInternalForce;
  const double    expectedB = forceOnMid / 4 + forceOnMidAB / 2 + unbalancedforceOnBC * 0.75 + 0.1 * unbalancedInternalForce;
  const double    expectedC = forceOnMid / 4 + forceOnC + unbalancedforceOnBC * 0.25 + 0.2 * unbalancedInternalForce;
  const double    expectedD = forceOnMid / 4 + 0.3 * unbalancedInternalForce;

  expected << expectedA, expectedB, expectedC, expectedD;
  // Check expected value, and conservation
  BOOST_CHECK(equals(expected, outValuesScalar));
  BOOST_CHECK(equals(netForce, outValuesScalar.sum()));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
