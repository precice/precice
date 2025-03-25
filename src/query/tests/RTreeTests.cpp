#include <Eigen/Core>
#include <algorithm>
#include <iterator>
#include <limits>
#include <list>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "logging/Logger.hpp"
#include "math/geometry.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "query/Index.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;
using namespace precice::query;

namespace {
PtrMesh fullMesh()
{
  PtrMesh ptr(new Mesh("MyMesh", 3, testing::nextMeshID()));
  auto   &mesh = *ptr;
  auto   &v1   = mesh.createVertex(Eigen::Vector3d(0, 2, 0));
  auto   &v2   = mesh.createVertex(Eigen::Vector3d(2, 1, 0));
  auto   &v3   = mesh.createVertex(Eigen::Vector3d(3, 0, 0));
  auto   &v4   = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  // Quad Borders
  auto &e1 = mesh.createEdge(v1, v2);
  auto &e2 = mesh.createEdge(v2, v3);
  auto &e3 = mesh.createEdge(v3, v4);
  auto &e4 = mesh.createEdge(v4, v1);
  // Diagonal
  auto &e5 = mesh.createEdge(v2, v4);
  // Triangles
  mesh.createTriangle(e1, e5, e4);
  mesh.createTriangle(e2, e3, e5);
  return ptr;
}

PtrMesh edgeMesh3D()
{
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 3, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  mesh->createVertex(Eigen::Vector3d(0, 0, 1));
  mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  mesh->createVertex(Eigen::Vector3d(0, 1, 1));
  mesh->createVertex(Eigen::Vector3d(1, 0, 0));
  mesh->createVertex(Eigen::Vector3d(1, 0, 1));
  auto &v1 = mesh->createVertex(Eigen::Vector3d(1, 1, 0));
  auto &v2 = mesh->createVertex(Eigen::Vector3d(1, 1, 1));
  mesh->createEdge(v1, v2);
  return mesh;
}

PtrMesh edgeMesh2D()
{
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector2d(0, 0));
  mesh->createVertex(Eigen::Vector2d(0, 1));
  auto &v1 = mesh->createVertex(Eigen::Vector2d(1, 0));
  auto &v2 = mesh->createVertex(Eigen::Vector2d(1, 1));
  mesh->createEdge(v1, v2);
  return mesh;
}

PtrMesh vertexMesh3D()
{
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 3, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  mesh->createVertex(Eigen::Vector3d(0, 0, 1));
  mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  mesh->createVertex(Eigen::Vector3d(0, 1, 1));
  mesh->createVertex(Eigen::Vector3d(1, 0, 0));
  mesh->createVertex(Eigen::Vector3d(1, 0, 1));
  mesh->createVertex(Eigen::Vector3d(1, 1, 0));
  mesh->createVertex(Eigen::Vector3d(1, 1, 1));
  mesh->computeBoundingBox();
  return mesh;
}
} // namespace

BOOST_AUTO_TEST_SUITE(QueryTests)
BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(Vertex)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Query2DVertex)
{
  PRECICE_TEST();
  auto            mesh = edgeMesh2D();
  Index           indexTree(mesh);
  Eigen::Vector2d location(0.2, 0.8);

  auto result = indexTree.getClosestVertex(location);
  BOOST_TEST(mesh->vertex(result.index).getCoords() == Eigen::Vector2d(0, 1));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Query3DVertex)
{
  PRECICE_TEST();
  auto            mesh = edgeMesh3D();
  Index           indexTree(mesh);
  Eigen::Vector3d location(0.8, 0.0, 0.8);

  auto result = indexTree.getClosestVertex(location);
  BOOST_TEST(mesh->vertex(result.index).getCoords() == Eigen::Vector3d(1, 0, 1));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Query3DFullVertex)
{
  PRECICE_TEST();
  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, precice::testing::nextMeshID()));
  const double z1  = 0.1;
  const double z2  = -0.1;
  auto        &v00 = mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  auto        &v01 = mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  auto        &v10 = mesh->createVertex(Eigen::Vector3d(1, 0, z1));
  auto        &v11 = mesh->createVertex(Eigen::Vector3d(1, 1, z1));
  auto        &v20 = mesh->createVertex(Eigen::Vector3d(2, 0, z2));
  auto        &v21 = mesh->createVertex(Eigen::Vector3d(2, 1, z2));
  auto        &v30 = mesh->createVertex(Eigen::Vector3d(3, 0, z2));
  auto        &v31 = mesh->createVertex(Eigen::Vector3d(3, 1, z2));
  auto        &ell = mesh->createEdge(v00, v01);
  auto        &elt = mesh->createEdge(v01, v11);
  auto        &elr = mesh->createEdge(v11, v10);
  auto        &elb = mesh->createEdge(v10, v00);
  auto        &eld = mesh->createEdge(v00, v11);
  auto        &erl = elr;
  auto        &ert = mesh->createEdge(v11, v21);
  auto        &err = mesh->createEdge(v21, v20);
  auto        &erb = mesh->createEdge(v20, v10);
  auto        &erd = mesh->createEdge(v10, v21);
  mesh->createTriangle(ell, elt, eld);
  mesh->createTriangle(eld, elb, elr);
  mesh->createTriangle(erl, ert, erd);
  mesh->createTriangle(erd, erb, err);

  Index indexTree(mesh);
  {
    Eigen::Vector3d location(0.8, 0.0, 0.8);
    auto            result = indexTree.getClosestVertex(location);

    BOOST_TEST(mesh->vertex(result.index).getID() == v10.getID());
  }
  {
    Eigen::Vector3d       location(0.8, 0.0, 0.8);
    int                   nVertices = 2;
    std::vector<VertexID> expectedResult({v00.getID(), v10.getID()});
    auto                  result = indexTree.getClosestVertices(location, nVertices);

    BOOST_TEST(result.size() == nVertices);
    BOOST_TEST(std::is_permutation(result.begin(), result.end(), expectedResult.begin()));
  }
  {
    Eigen::Vector3d       location(3.5, 3.5, 0.0);
    int                   nVertices = 4;
    std::vector<VertexID> expectedResult({v11.getID(), v30.getID(), v21.getID(), v31.getID()});
    auto                  result = indexTree.getClosestVertices(location, nVertices);

    BOOST_TEST(result.size() == nVertices);
    BOOST_TEST(std::is_permutation(result.begin(), result.end(), expectedResult.begin()));
  }
}

/// Resembles how boost geometry is used inside the PetRBF
PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(QueryWithBoxEmpty)
{
  PRECICE_TEST();
  auto         mesh = vertexMesh3D();
  Index        indexTree(mesh);
  mesh::Vertex searchVertex(Eigen::Vector3d(0.8, 1, 0), 0);
  double       radius = 0.1; // No vertices in radius

  auto results = indexTree.getVerticesInsideBox(searchVertex, radius);
  BOOST_TEST(results.empty());
}

/// Resembles how boost geometry is used inside the PetRBF
PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(QueryWithBox2Matches)
{
  PRECICE_TEST();
  auto  mesh = vertexMesh3D();
  Index indexTree(mesh);

  mesh::Vertex searchVertex(Eigen::Vector3d(0.8, 1, 0), 0);
  double       radius = 0.81; // Two vertices in radius

  auto results = indexTree.getVerticesInsideBox(searchVertex, radius);
  BOOST_TEST(results.size() == 2);
  BOOST_TEST(mesh->vertex(results.at(0)).getCoords() == Eigen::Vector3d(0, 1, 0));
  BOOST_TEST(mesh->vertex(results.at(1)).getCoords() == Eigen::Vector3d(1, 1, 0));
}

/// Resembles how boost geometry is used inside the PetRBF
PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(QueryWithBoxEverything)
{
  PRECICE_TEST();
  auto  mesh = vertexMesh3D();
  Index indexTree(mesh);

  mesh::Vertex searchVertex(Eigen::Vector3d(0.8, 1, 0), 0);
  double       radius = std::numeric_limits<double>::max();

  auto results = indexTree.getVerticesInsideBox(searchVertex, radius);
  BOOST_TEST(results.size() == 8);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(QueryRtreeBoundingBox2D)
{
  PRECICE_TEST();
  auto mesh   = edgeMesh2D();
  auto result = mesh->index().getRtreeBounds();
  mesh->computeBoundingBox();
  auto comparison = mesh->getBoundingBox();

  BOOST_TEST(result == comparison);
  BOOST_TEST(result.minCorner() == Eigen::Vector2d(0, 0));
  BOOST_TEST(result.maxCorner() == Eigen::Vector2d(1, 1));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(QueryRtreeBoundingBox3D)
{
  PRECICE_TEST();
  auto mesh   = vertexMesh3D();
  auto result = mesh->index().getRtreeBounds();
  mesh->computeBoundingBox();
  auto comparison = mesh->getBoundingBox();

  BOOST_TEST(result == comparison);
  BOOST_TEST(result.minCorner() == Eigen::Vector3d(0, 0, 0));
  BOOST_TEST(result.maxCorner() == Eigen::Vector3d(1, 1, 1));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(QueryRtreeBoundingBox3DComplex)
{
  PRECICE_TEST();
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 3, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector3d(7, 4, 3.3));
  mesh->createVertex(Eigen::Vector3d(26.4777, 5, 8));
  mesh->createVertex(Eigen::Vector3d(-23.4, 100000.2, 7));
  mesh->createVertex(Eigen::Vector3d(0.211, -21.37, 0.00003));
  auto result = mesh->index().getRtreeBounds();
  mesh->computeBoundingBox();
  auto comparison = mesh->getBoundingBox();

  BOOST_TEST(result == comparison);
  BOOST_TEST(result.minCorner() == Eigen::Vector3d(-23.4, -21.37, 0.00003));
  BOOST_TEST(result.maxCorner() == Eigen::Vector3d(26.4777, 100000.2, 8));
}

BOOST_AUTO_TEST_SUITE_END() // Vertex

BOOST_AUTO_TEST_SUITE(Edge)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Query2DEdge)
{
  PRECICE_TEST();
  auto            mesh = edgeMesh2D();
  Index           indexTree(mesh);
  Eigen::Vector2d location(0.2, 0.8);

  auto results = indexTree.getClosestEdges(location, 1);
  BOOST_TEST(results.size() == 1);
  auto &edge = mesh->edges().at(results.front().index);

  BOOST_TEST(edge.vertex(0).getCoords() == Eigen::Vector2d(1, 0));
  BOOST_TEST(edge.vertex(1).getCoords() == Eigen::Vector2d(1, 1));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Query3DEdge)
{
  PRECICE_TEST();
  auto            mesh = edgeMesh3D();
  Index           indexTree(mesh);
  Eigen::Vector3d location(1.8, 0.0, 0.8);

  auto results = indexTree.getClosestEdges(location, 1);

  BOOST_TEST(results.size() == 1);
  auto match = results.front().index;

  BOOST_TEST(match < mesh->edges().size());
  auto           &edge = mesh->edges().at(match);
  Eigen::Vector3d p1(1, 1, 0);
  Eigen::Vector3d p2(1, 1, 1);
  BOOST_TEST((edge.vertex(0).getCoords() == p1 || edge.vertex(0).getCoords() == p2));
  if (edge.vertex(0).getCoords() == p1) {
    BOOST_TEST(edge.vertex(1).getCoords() == p2);
  } else {
    BOOST_TEST(edge.vertex(1).getCoords() == p1);
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Query3DFullEdge)
{
  PRECICE_TEST();
  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, precice::testing::nextMeshID()));
  const double z1  = 0.1;
  const double z2  = -0.1;
  auto        &v00 = mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  auto        &v01 = mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  auto        &v10 = mesh->createVertex(Eigen::Vector3d(1, 0, z1));
  auto        &v11 = mesh->createVertex(Eigen::Vector3d(1, 1, z1));
  auto        &v20 = mesh->createVertex(Eigen::Vector3d(2, 0, z2));
  auto        &v21 = mesh->createVertex(Eigen::Vector3d(2, 1, z2));
  auto        &ell = mesh->createEdge(v00, v01);
  auto        &elt = mesh->createEdge(v01, v11);
  auto        &elr = mesh->createEdge(v11, v10);
  auto        &elb = mesh->createEdge(v10, v00);
  auto        &eld = mesh->createEdge(v00, v11);
  auto        &erl = elr;
  auto        &ert = mesh->createEdge(v11, v21);
  auto        &err = mesh->createEdge(v21, v20);
  auto        &erb = mesh->createEdge(v20, v10);
  auto        &erd = mesh->createEdge(v10, v21);
  mesh->createTriangle(ell, elt, eld);
  mesh->createTriangle(eld, elb, elr);
  mesh->createTriangle(erl, ert, erd);
  mesh->createTriangle(erd, erb, err);

  Index           indexTree(mesh);
  Eigen::Vector3d location(0.8, 0.5, 0.0);
  auto            results = indexTree.getClosestEdges(location, 2);

  std::set<mesh::Edge *> matches{
      &mesh->edges().at(results.at(0).index),
      &mesh->edges().at(results.at(1).index),
  };
  std::set<mesh::Edge *> expected{&elr, &eld};

  BOOST_TEST(matches.size() == 2);
  BOOST_TEST(matches == expected, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_SUITE_END() // Edge

BOOST_AUTO_TEST_SUITE(Triangle)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Query3DFullTriangle)
{
  PRECICE_TEST();

  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, precice::testing::nextMeshID()));
  const double z1  = 0.1;
  const double z2  = -0.1;
  auto        &v00 = mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  auto        &v01 = mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  auto        &v10 = mesh->createVertex(Eigen::Vector3d(1, 0, z1));
  auto        &v11 = mesh->createVertex(Eigen::Vector3d(1, 1, z1));
  auto        &v20 = mesh->createVertex(Eigen::Vector3d(2, 0, z2));
  auto        &v21 = mesh->createVertex(Eigen::Vector3d(2, 1, z2));
  auto        &ell = mesh->createEdge(v00, v01);
  auto        &elt = mesh->createEdge(v01, v11);
  auto        &elr = mesh->createEdge(v11, v10);
  auto        &elb = mesh->createEdge(v10, v00);
  auto        &eld = mesh->createEdge(v00, v11);
  auto        &erl = elr;
  auto        &ert = mesh->createEdge(v11, v21);
  auto        &err = mesh->createEdge(v21, v20);
  auto        &erb = mesh->createEdge(v20, v10);
  auto        &erd = mesh->createEdge(v10, v21);
  auto        &tlt = mesh->createTriangle(ell, elt, eld);
  auto        &tlb = mesh->createTriangle(eld, elb, elr);
  auto        &trt = mesh->createTriangle(erl, ert, erd);
  mesh->createTriangle(erd, erb, err);

  Index indexTree(mesh);

  Eigen::Vector3d location(0.7, 0.5, 0.0);

  auto results = indexTree.getClosestTriangles(location, 3);
  BOOST_TEST(results.size() == 3);

  std::set<mesh::Triangle *> matches{
      &mesh->triangles().at(results.at(0).index),
      &mesh->triangles().at(results.at(1).index),
      &mesh->triangles().at(results.at(2).index)};
  std::set<mesh::Triangle *> expected{&tlb, &tlt, &trt};

  BOOST_TEST(matches.size() == 3);
  BOOST_TEST(matches == expected, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_SUITE_END() // Triangle

BOOST_AUTO_TEST_SUITE(Projection)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ProjectionToVertex)
{
  PRECICE_TEST();
  auto  meshPtr = fullMesh();
  Index indexTree(meshPtr);

  Eigen::Vector3d     location(4.0, 0.0, 0.0);
  std::vector<int>    expectedIndices = {2};
  std::vector<double> expectedWeights = {1.0};

  auto match = indexTree.findNearestProjection(location, 1);

  BOOST_TEST(match.polation.getWeightedElements().size() == 1); // Check number of weights
  BOOST_TEST(match.polation.distance() == 1.0);                 // Check the distance
  BOOST_TEST(match.polation.isInterpolation());

  for (int i = 0; i < static_cast<int>(match.polation.getWeightedElements().size()); ++i) {
    BOOST_TEST(match.polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i)); // Check index
    BOOST_TEST(match.polation.getWeightedElements().at(i).weight == expectedWeights.at(i));   // Check the weight
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ProjectionToEdge)
{
  PRECICE_TEST();
  auto  meshPtr = fullMesh();
  Index indexTree(meshPtr);

  Eigen::Vector3d     location(2.0, -1.0, 0.0);
  std::vector<int>    expectedIndices = {2, 3};
  std::vector<double> expectedWeights = {0.5, 0.5};

  auto match = indexTree.findNearestProjection(location, 1);

  BOOST_TEST(match.polation.getWeightedElements().size() == 2); // Check number of weights
  BOOST_TEST(match.polation.distance() == 1.0);                 // Check the distance
  BOOST_TEST(match.polation.isInterpolation());

  for (int i = 0; i < static_cast<int>(match.polation.getWeightedElements().size()); ++i) {
    BOOST_TEST(match.polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i)); // Check index
    BOOST_TEST(match.polation.getWeightedElements().at(i).weight == expectedWeights.at(i));   // Check the weight
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ProjectionToTriangle)
{
  PRECICE_TEST();
  auto  meshPtr = fullMesh();
  Index indexTree(meshPtr);

  Eigen::Vector3d     location(1.0, 1.0, 0.1);
  std::vector<int>    expectedIndices = {0, 1, 3};
  std::vector<double> expectedWeights = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

  auto match = indexTree.findNearestProjection(location, 1);

  BOOST_TEST(match.polation.getWeightedElements().size() == 3); // Check number of weights
  BOOST_TEST(match.polation.distance() == 0.1);                 // Check the distance
  BOOST_TEST(match.polation.isInterpolation());

  for (int i = 0; i < static_cast<int>(match.polation.getWeightedElements().size()); ++i) {
    BOOST_TEST(match.polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i)); // Check index
    BOOST_TEST(match.polation.getWeightedElements().at(i).weight == expectedWeights.at(i));   // Check the weight
  }
}

BOOST_AUTO_TEST_SUITE_END() // Projection

BOOST_AUTO_TEST_SUITE(Tetrahedra)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(CubeBoundingBoxIndex)
{
  PRECICE_TEST();
  PtrMesh ptr(new Mesh("MyMesh", 3, testing::nextMeshID()));
  auto   &mesh = *ptr;
  Index   indexTree(ptr);

  Eigen::Vector3d  location(0.5, 0.5, 0.5);
  std::vector<int> expectedIndices = {0, 1};
  // Set up 2 tetrahedra with the same bounding box
  auto &v00 = mesh.createVertex(Eigen::Vector3d(0, 0, 0));
  auto &v01 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto &v02 = mesh.createVertex(Eigen::Vector3d(0, 1, 0));
  auto &v03 = mesh.createVertex(Eigen::Vector3d(1, 0, 1));
  auto &v04 = mesh.createVertex(Eigen::Vector3d(1, 1, 1));

  mesh.createTetrahedron(v00, v01, v02, v03);
  mesh.createTetrahedron(v04, v01, v02, v03);

  auto match = indexTree.getEnclosingTetrahedra(location);

  BOOST_TEST(match.size() == 2);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TetraIndexing)
{
  PRECICE_TEST();
  /*
  For a location and 3 tetrahedra such that:
  - First contains the location
  - Second doesn't, but its Bounding Box does
  - Third doesn't and neither does its AABB
  Check that only 1st and 2nd are found by getEnclosingTetrahedra
  */
  PtrMesh ptr(new Mesh("MyMesh", 3, testing::nextMeshID()));
  auto   &mesh = *ptr;
  Index   indexTree(ptr);

  Eigen::Vector3d  location(0.2, 0.2, 0.2);
  std::vector<int> expectedIndices = {0, 1};

  // Set containing tetra
  auto &v00 = mesh.createVertex(Eigen::Vector3d(0, 0, 0));
  auto &v01 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto &v02 = mesh.createVertex(Eigen::Vector3d(0, 1, 0));
  auto &v03 = mesh.createVertex(Eigen::Vector3d(0, 0, 1));
  mesh.createTetrahedron(v00, v01, v02, v03);

  // Set non-containing tetra with containing BB (from 0 to 1 in each direction)
  auto &v10 = mesh.createVertex(Eigen::Vector3d(1, 1, 1));
  auto &v11 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto &v12 = mesh.createVertex(Eigen::Vector3d(0, 1, 0));
  auto &v13 = mesh.createVertex(Eigen::Vector3d(0, 0, 1));
  mesh.createTetrahedron(v10, v11, v12, v13);

  // Set tetra far away
  auto &v20 = mesh.createVertex(Eigen::Vector3d(1, 1, 1));
  auto &v21 = mesh.createVertex(Eigen::Vector3d(2, 1, 1));
  auto &v22 = mesh.createVertex(Eigen::Vector3d(1, 2, 1));
  auto &v23 = mesh.createVertex(Eigen::Vector3d(1, 1, 2));
  mesh.createTetrahedron(v20, v21, v22, v23);

  auto match = indexTree.getEnclosingTetrahedra(location);

  BOOST_TEST(match.size() == 2);
  BOOST_TEST(((match[0] == 0 && match[1] == 1) || (match[0] == 1 && match[1] == 0)));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TetraWorksOnBoundary)
{
  PRECICE_TEST();
  /*
 Check that the AABB safety factor is high enough. Do all the corners of a tetra fit inside its AABB?
  */
  PtrMesh ptr(new Mesh("MyMesh", 3, testing::nextMeshID()));
  auto   &mesh = *ptr;
  Index   indexTree(ptr);

  std::vector<Eigen::Vector3d> locations;
  locations.emplace_back(0, 0, 0);
  locations.emplace_back(1, 0, 0);
  locations.emplace_back(0, 1, 0);
  locations.emplace_back(0, 0, 1);

  // Set containing tetra
  auto &v00 = mesh.createVertex(Eigen::Vector3d(0, 0, 0));
  auto &v01 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto &v02 = mesh.createVertex(Eigen::Vector3d(0, 1, 0));
  auto &v03 = mesh.createVertex(Eigen::Vector3d(0, 0, 1));
  mesh.createTetrahedron(v00, v01, v02, v03);

  for (const auto &vertex : locations) {
    auto match = indexTree.getEnclosingTetrahedra(vertex);
    BOOST_TEST(match.size() == 1);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Tetrahedra

BOOST_AUTO_TEST_SUITE_END() // Mesh
BOOST_AUTO_TEST_SUITE_END() // Query
