#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "logging/LogMacros.hpp"
#include "mapping/LinearCellInterpolation.hpp"
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
BOOST_AUTO_TEST_SUITE(LinearCellInterpolation)

BOOST_AUTO_TEST_CASE(ConsistentNonIncremental)
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
  Should pass: triangles are discarded when the mapping doesn't require the in SolverInterface::setMeshTriangle,
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
  outMesh->createVertex(Eigen::Vector2d(-0.1, -0.1)); // Currently maps to opposite edge
  outMesh->createVertex(Eigen::Vector2d(2.5, -1.0));
  outMesh->createVertex(Eigen::Vector2d(2.5, 10.0));

  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::LinearCellInterpolation mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check expected
  Eigen::VectorXd expected(outMesh->vertices().size());
  expected << 2.0, 1.0, 2.0, 3.0, 1.49, 2.25, 1.5, 2.5, 2.0, 3.0;
  BOOST_CHECK(equals(expected, outValuesScalar));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
