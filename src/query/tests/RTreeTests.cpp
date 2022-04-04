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
#include "query/impl/Indexer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;
using namespace precice::query;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace {
PtrMesh fullMesh()
{
  PtrMesh ptr(new Mesh("MyMesh", 3, testing::nextMeshID()));
  auto &  mesh = *ptr;
  auto &  v1   = mesh.createVertex(Eigen::Vector3d(0, 2, 0));
  auto &  v2   = mesh.createVertex(Eigen::Vector3d(2, 1, 0));
  auto &  v3   = mesh.createVertex(Eigen::Vector3d(3, 0, 0));
  auto &  v4   = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
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

BOOST_AUTO_TEST_CASE(Query2DVertex)
{
  PRECICE_TEST(1_rank);
  auto            mesh = edgeMesh2D();
  Index           indexTree(mesh);
  Eigen::Vector2d location(0.2, 0.8);

  auto result = indexTree.getClosestVertex(location);
  BOOST_TEST(mesh->vertices().at(result.index).getCoords() == Eigen::Vector2d(0, 1));
  BOOST_TEST(result.distance == 0.28284271247461906);
}

BOOST_AUTO_TEST_CASE(Query3DVertex)
{
  PRECICE_TEST(1_rank);
  auto            mesh = edgeMesh3D();
  Index           indexTree(mesh);
  Eigen::Vector3d location(0.8, 0.0, 0.8);

  auto result = indexTree.getClosestVertex(location);
  BOOST_TEST(mesh->vertices().at(result.index).getCoords() == Eigen::Vector3d(1, 0, 1));
  BOOST_TEST(result.distance == 0.28284271247461906);
}

BOOST_AUTO_TEST_CASE(Query3DFullVertex)
{
  PRECICE_TEST(1_rank);
  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, precice::testing::nextMeshID()));
  const double z1  = 0.1;
  const double z2  = -0.1;
  auto &       v00 = mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  auto &       v01 = mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  auto &       v10 = mesh->createVertex(Eigen::Vector3d(1, 0, z1));
  auto &       v11 = mesh->createVertex(Eigen::Vector3d(1, 1, z1));
  auto &       v20 = mesh->createVertex(Eigen::Vector3d(2, 0, z2));
  auto &       v21 = mesh->createVertex(Eigen::Vector3d(2, 1, z2));
  auto &       ell = mesh->createEdge(v00, v01);
  auto &       elt = mesh->createEdge(v01, v11);
  auto &       elr = mesh->createEdge(v11, v10);
  auto &       elb = mesh->createEdge(v10, v00);
  auto &       eld = mesh->createEdge(v00, v11);
  auto &       erl = elr;
  auto &       ert = mesh->createEdge(v11, v21);
  auto &       err = mesh->createEdge(v21, v20);
  auto &       erb = mesh->createEdge(v20, v10);
  auto &       erd = mesh->createEdge(v10, v21);
  mesh->createTriangle(ell, elt, eld);
  mesh->createTriangle(eld, elb, elr);
  mesh->createTriangle(erl, ert, erd);
  mesh->createTriangle(erd, erb, err);

  Eigen::Vector3d location(0.8, 0.0, 0.8);
  Index           indexTree(mesh);
  auto            result = indexTree.getClosestVertex(location);

  BOOST_TEST(mesh->vertices().at(result.index).getID() == v10.getID());
}

/// Resembles how boost geometry is used inside the PetRBF
BOOST_AUTO_TEST_CASE(QueryWithBoxEmpty)
{
  PRECICE_TEST(1_rank);
  auto         mesh = vertexMesh3D();
  Index        indexTree(mesh);
  mesh::Vertex searchVertex(Eigen::Vector3d(0.8, 1, 0), 0);
  double       radius = 0.1; // No vertices in radius

  auto results = indexTree.getVerticesInsideBox(searchVertex, radius);
  BOOST_TEST(results.empty());
}

/// Resembles how boost geometry is used inside the PetRBF
BOOST_AUTO_TEST_CASE(QueryWithBox2Matches)
{
  PRECICE_TEST(1_rank);
  auto  mesh = vertexMesh3D();
  Index indexTree(mesh);
  auto  tree = impl::Indexer::instance()->getVertexRTree(mesh);
  BOOST_TEST(tree->size() == 8);

  mesh::Vertex searchVertex(Eigen::Vector3d(0.8, 1, 0), 0);
  double       radius = 0.81; // Two vertices in radius

  auto results = indexTree.getVerticesInsideBox(searchVertex, radius);
  BOOST_TEST(results.size() == 2);
  BOOST_TEST(mesh->vertices().at(results.at(0)).getCoords() == Eigen::Vector3d(0, 1, 0));
  BOOST_TEST(mesh->vertices().at(results.at(1)).getCoords() == Eigen::Vector3d(1, 1, 0));
}

/// Resembles how boost geometry is used inside the PetRBF
BOOST_AUTO_TEST_CASE(QueryWithBoxEverything)
{
  PRECICE_TEST(1_rank);
  auto  mesh = vertexMesh3D();
  Index indexTree(mesh);
  auto  tree = impl::Indexer::instance()->getVertexRTree(mesh);
  BOOST_TEST(tree->size() == 8);

  mesh::Vertex searchVertex(Eigen::Vector3d(0.8, 1, 0), 0);
  double       radius = std::numeric_limits<double>::max();

  auto results = indexTree.getVerticesInsideBox(searchVertex, radius);
  BOOST_TEST(results.size() == 8);
}

BOOST_AUTO_TEST_SUITE_END() // Vertex

BOOST_AUTO_TEST_SUITE(Edge)

BOOST_AUTO_TEST_CASE(Query2DEdge)
{
  PRECICE_TEST(1_rank);
  auto            mesh = edgeMesh2D();
  Index           indexTree(mesh);
  Eigen::Vector2d location(0.2, 0.8);

  auto results = indexTree.getClosestEdges(location, 1);
  BOOST_TEST(results.size() == 1);
  auto &edge = mesh->edges().at(results.front().index);

  BOOST_TEST(edge.vertex(0).getCoords() == Eigen::Vector2d(1, 0));
  BOOST_TEST(edge.vertex(1).getCoords() == Eigen::Vector2d(1, 1));
}

BOOST_AUTO_TEST_CASE(Query3DEdge)
{
  PRECICE_TEST(1_rank);
  auto            mesh = edgeMesh3D();
  Index           indexTree(mesh);
  Eigen::Vector3d location(1.8, 0.0, 0.8);

  auto results = indexTree.getClosestEdges(location, 1);

  BOOST_TEST(results.size() == 1);
  auto match = results.front().index;

  BOOST_TEST(match < mesh->edges().size());
  auto &          edge = mesh->edges().at(match);
  Eigen::Vector3d p1(1, 1, 0);
  Eigen::Vector3d p2(1, 1, 1);
  BOOST_TEST((edge.vertex(0).getCoords() == p1 || edge.vertex(0).getCoords() == p2));
  if (edge.vertex(0).getCoords() == p1) {
    BOOST_TEST(edge.vertex(1).getCoords() == p2);
  } else {
    BOOST_TEST(edge.vertex(1).getCoords() == p1);
  }
}

BOOST_AUTO_TEST_CASE(Query3DFullEdge)
{
  PRECICE_TEST(1_rank);
  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, precice::testing::nextMeshID()));
  const double z1  = 0.1;
  const double z2  = -0.1;
  auto &       v00 = mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  auto &       v01 = mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  auto &       v10 = mesh->createVertex(Eigen::Vector3d(1, 0, z1));
  auto &       v11 = mesh->createVertex(Eigen::Vector3d(1, 1, z1));
  auto &       v20 = mesh->createVertex(Eigen::Vector3d(2, 0, z2));
  auto &       v21 = mesh->createVertex(Eigen::Vector3d(2, 1, z2));
  auto &       ell = mesh->createEdge(v00, v01);
  auto &       elt = mesh->createEdge(v01, v11);
  auto &       elr = mesh->createEdge(v11, v10);
  auto &       elb = mesh->createEdge(v10, v00);
  auto &       eld = mesh->createEdge(v00, v11);
  auto &       erl = elr;
  auto &       ert = mesh->createEdge(v11, v21);
  auto &       err = mesh->createEdge(v21, v20);
  auto &       erb = mesh->createEdge(v20, v10);
  auto &       erd = mesh->createEdge(v10, v21);
  mesh->createTriangle(ell, elt, eld);
  mesh->createTriangle(eld, elb, elr);
  mesh->createTriangle(erl, ert, erd);
  mesh->createTriangle(erd, erb, err);

  Index           indexTree(mesh);
  Eigen::Vector3d location(0.8, 0.5, 0.0);
  auto            results = indexTree.getClosestEdges(location, 2);

  BOOST_TEST(results.size() == 2);
  BOOST_TEST(results.at(0).index == eld.getID());
  BOOST_TEST(results.at(1).index == elr.getID());
}

BOOST_AUTO_TEST_SUITE_END() // Edge

BOOST_AUTO_TEST_SUITE(Triangle)

BOOST_AUTO_TEST_CASE(Query3DFullTriangle)
{
  PRECICE_TEST(1_rank);

  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, precice::testing::nextMeshID()));
  const double z1  = 0.1;
  const double z2  = -0.1;
  auto &       v00 = mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  auto &       v01 = mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  auto &       v10 = mesh->createVertex(Eigen::Vector3d(1, 0, z1));
  auto &       v11 = mesh->createVertex(Eigen::Vector3d(1, 1, z1));
  auto &       v20 = mesh->createVertex(Eigen::Vector3d(2, 0, z2));
  auto &       v21 = mesh->createVertex(Eigen::Vector3d(2, 1, z2));
  auto &       ell = mesh->createEdge(v00, v01);
  auto &       elt = mesh->createEdge(v01, v11);
  auto &       elr = mesh->createEdge(v11, v10);
  auto &       elb = mesh->createEdge(v10, v00);
  auto &       eld = mesh->createEdge(v00, v11);
  auto &       erl = elr;
  auto &       ert = mesh->createEdge(v11, v21);
  auto &       err = mesh->createEdge(v21, v20);
  auto &       erb = mesh->createEdge(v20, v10);
  auto &       erd = mesh->createEdge(v10, v21);
  auto &       tlt = mesh->createTriangle(ell, elt, eld);
  auto &       tlb = mesh->createTriangle(eld, elb, elr);
  auto &       trt = mesh->createTriangle(erl, ert, erd);
  auto &       trb = mesh->createTriangle(erd, erb, err);

  Index indexTree(mesh);

  Eigen::Vector3d location(0.7, 0.5, 0.0);

  auto results = indexTree.getClosestTriangles(location, 3);
  BOOST_TEST(results.size() == 3);
  BOOST_TEST(results.at(0).index == tlb.getID());
  BOOST_TEST(results.at(1).index == tlt.getID());
  BOOST_TEST(results.at(2).index == trt.getID());
  BOOST_TEST(results.at(2).index != trb.getID());
}

BOOST_AUTO_TEST_SUITE_END() // Triangle

BOOST_AUTO_TEST_SUITE(Cache)

BOOST_AUTO_TEST_CASE(ClearOnChange)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector2d(0, 0));

  // The Cache should clear whenever a mesh changes
  auto vTree = query::impl::Indexer::instance()->getVertexRTree(mesh);
  BOOST_TEST(query::impl::Indexer::instance()->getCacheSize() == 1);
  mesh->meshChanged(*mesh); // Emit signal, that mesh has changed
  BOOST_TEST(query::impl::Indexer::instance()->getCacheSize() == 0);
}

BOOST_AUTO_TEST_CASE(ClearOnDestruction)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector2d(0, 0));

  // The Cache should clear whenever we destroy the Mesh
  auto vTree = query::impl::Indexer::instance()->getVertexRTree(mesh);
  BOOST_TEST(query::impl::Indexer::instance()->getCacheSize() == 1);
  mesh.reset(); // Destroy mesh object, signal is emitted to clear cache
  BOOST_TEST(query::impl::Indexer::instance()->getCacheSize() == 0);
}

BOOST_AUTO_TEST_CASE(CacheVertices)
{
  PRECICE_TEST(1_rank);
  auto ptr = fullMesh();

  auto vt1 = impl::Indexer::instance()->getVertexRTree(ptr);
  auto vt2 = impl::Indexer::instance()->getVertexRTree(ptr);
  BOOST_TEST(vt1 == vt2);
}

BOOST_AUTO_TEST_CASE(CacheEdges)
{
  PRECICE_TEST(1_rank);
  auto ptr = fullMesh();

  auto et1 = impl::Indexer::instance()->getEdgeRTree(ptr);
  auto et2 = impl::Indexer::instance()->getEdgeRTree(ptr);
  BOOST_TEST(et1 == et2);
}

BOOST_AUTO_TEST_CASE(CacheTriangles)
{
  PRECICE_TEST(1_rank);
  auto ptr = fullMesh();

  auto tt1 = impl::Indexer::instance()->getTriangleRTree(ptr);
  auto tt2 = impl::Indexer::instance()->getTriangleRTree(ptr);
  BOOST_TEST(tt1 == tt2);
}

BOOST_AUTO_TEST_CASE(CacheAll)
{
  PRECICE_TEST(1_rank);
  auto ptr = fullMesh();

  auto vt1 = impl::Indexer::instance()->getVertexRTree(ptr);
  auto et1 = impl::Indexer::instance()->getEdgeRTree(ptr);
  auto tt1 = impl::Indexer::instance()->getTriangleRTree(ptr);

  auto vt2 = impl::Indexer::instance()->getVertexRTree(ptr);
  auto et2 = impl::Indexer::instance()->getEdgeRTree(ptr);
  auto tt2 = impl::Indexer::instance()->getTriangleRTree(ptr);

  BOOST_TEST(vt1 == vt2);
  BOOST_TEST(et1 == et2);
  BOOST_TEST(tt1 == tt2);
}

BOOST_AUTO_TEST_SUITE_END() // Cache

BOOST_AUTO_TEST_SUITE(Projection)

BOOST_AUTO_TEST_CASE(ProjectionToVertex)
{
  PRECICE_TEST(1_rank);
  auto  meshPtr = fullMesh();
  Index indexTree(meshPtr);

  Eigen::Vector3d     location(4.0, 0.0, 0.0);
  std::vector<int>    expectedIndices = {2};
  std::vector<double> expectedWeights = {1.0};

  auto match = indexTree.findNearestProjection(location, 1);

  BOOST_TEST(match.polation.getWeightedElements().size() == 1); // Check number of weights
  BOOST_TEST(match.distance == 1.0);                            // Check the distance
  BOOST_TEST(match.polation.isInterpolation());

  for (int i = 0; i < static_cast<int>(match.polation.getWeightedElements().size()); ++i) {
    BOOST_TEST(match.polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i)); // Check index
    BOOST_TEST(match.polation.getWeightedElements().at(i).weight == expectedWeights.at(i));   // Check the weight
  }
}

BOOST_AUTO_TEST_CASE(ProjectionToEdge)
{
  PRECICE_TEST(1_rank);
  auto  meshPtr = fullMesh();
  Index indexTree(meshPtr);

  Eigen::Vector3d     location(2.0, -1.0, 0.0);
  std::vector<int>    expectedIndices = {2, 3};
  std::vector<double> expectedWeights = {0.5, 0.5};

  auto match = indexTree.findNearestProjection(location, 1);

  BOOST_TEST(match.polation.getWeightedElements().size() == 2); // Check number of weights
  BOOST_TEST(match.distance == 1.0);                            // Check the distance
  BOOST_TEST(match.polation.isInterpolation());

  for (int i = 0; i < static_cast<int>(match.polation.getWeightedElements().size()); ++i) {
    BOOST_TEST(match.polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i)); // Check index
    BOOST_TEST(match.polation.getWeightedElements().at(i).weight == expectedWeights.at(i));   // Check the weight
  }
}

BOOST_AUTO_TEST_CASE(ProjectionToTriangle)
{
  PRECICE_TEST(1_rank);
  auto  meshPtr = fullMesh();
  Index indexTree(meshPtr);

  Eigen::Vector3d     location(1.0, 1.0, 0.1);
  std::vector<int>    expectedIndices = {0, 1, 3};
  std::vector<double> expectedWeights = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

  auto match = indexTree.findNearestProjection(location, 1);

  BOOST_TEST(match.polation.getWeightedElements().size() == 3); // Check number of weights
  BOOST_TEST(match.distance == 0.0);                            // Check the distance
  BOOST_TEST(match.polation.isInterpolation());

  for (int i = 0; i < static_cast<int>(match.polation.getWeightedElements().size()); ++i) {
    BOOST_TEST(match.polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i)); // Check index
    BOOST_TEST(match.polation.getWeightedElements().at(i).weight == expectedWeights.at(i));   // Check the weight
  }
}

BOOST_AUTO_TEST_SUITE_END() // Projection

BOOST_AUTO_TEST_SUITE_END() // Mesh
BOOST_AUTO_TEST_SUITE_END() // Query
