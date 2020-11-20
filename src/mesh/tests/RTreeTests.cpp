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
#include "mesh/RTree.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/impl/BBUtils.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace {
PtrMesh fullMesh()
{
  PtrMesh ptr(new Mesh("MyMesh", 3, false, testing::nextMeshID()));
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
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 3, false, precice::testing::nextMeshID()));
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
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, false, testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector2d(0, 0));
  mesh->createVertex(Eigen::Vector2d(0, 1));
  auto &v1 = mesh->createVertex(Eigen::Vector2d(1, 0));
  auto &v2 = mesh->createVertex(Eigen::Vector2d(1, 1));
  mesh->createEdge(v1, v2);
  return mesh;
}

PtrMesh vertexMesh3D()
{
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 3, false, precice::testing::nextMeshID()));
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

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(RTree)

BOOST_AUTO_TEST_SUITE(Query)

BOOST_AUTO_TEST_SUITE(Vertex)

BOOST_AUTO_TEST_CASE(Query2DVertex)
{
  PRECICE_TEST(1_rank);
  auto mesh = edgeMesh2D();
  auto tree = rtree::getVertexRTree(mesh);

  BOOST_TEST(tree->size() == 4);

  Eigen::VectorXd     searchVector(Eigen::Vector2d(0.2, 0.8));
  std::vector<size_t> results;

  tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

  BOOST_TEST(results.size() == 1);
  BOOST_TEST(mesh->vertices().at(results.at(0)).getCoords() == Eigen::Vector2d(0, 1));
}

BOOST_AUTO_TEST_CASE(Query3DVertex)
{
  PRECICE_TEST(1_rank);
  auto mesh = edgeMesh3D();
  {
    auto tree = rtree::getVertexRTree(mesh);

    BOOST_TEST(tree->size() == 8);

    Eigen::VectorXd     searchVector(Eigen::Vector3d(0.8, 0.0, 0.8));
    std::vector<size_t> results;

    tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

    BOOST_TEST(results.size() == 1);
    BOOST_TEST(mesh->vertices().at(results.at(0)).getCoords() == Eigen::Vector3d(1, 0, 1));
  }
}

BOOST_AUTO_TEST_CASE(Query3DFullVertex)
{
  PRECICE_TEST(1_rank);
  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, false, precice::testing::nextMeshID()));
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

  auto tree = rtree::getVertexRTree(mesh);

  BOOST_TEST(tree->size() == 6);

  Eigen::VectorXd     searchVector(Eigen::Vector3d(0.8, 0.0, 0.8));
  std::vector<size_t> results;

  tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

  BOOST_TEST_INFO(results);
  BOOST_TEST(results.size() == 1);
  BOOST_TEST(mesh->vertices().at(results.front()).getID() == v10.getID());
}

/// Resembles how boost geometry is used inside the PetRBF
BOOST_AUTO_TEST_CASE(QueryWithBoxEmpty)
{
  PRECICE_TEST(1_rank);
  auto mesh = vertexMesh3D();

  auto tree = rtree::getVertexRTree(mesh);
  BOOST_TEST(tree->size() == 8);

  Eigen::VectorXd searchVector(Eigen::Vector3d(0.8, 1, 0));

  std::vector<size_t>             results;
  double                          radius = 0.1; // No vertices in radius
  bg::model::box<Eigen::VectorXd> search_box(
      searchVector - Eigen::VectorXd::Constant(3, radius),
      searchVector + Eigen::VectorXd::Constant(3, radius));

  tree->query(bg::index::intersects(search_box) and bg::index::satisfies([&](size_t const i) { return bg::distance(searchVector, mesh->vertices().at(i)) <= radius; }),
              std::back_inserter(results));

  BOOST_TEST(results.empty());
}

/// Resembles how boost geometry is used inside the PetRBF
BOOST_AUTO_TEST_CASE(QueryWithBox2Matches)
{
  PRECICE_TEST(1_rank);
  auto mesh = vertexMesh3D();

  auto tree = rtree::getVertexRTree(mesh);
  BOOST_TEST(tree->size() == 8);

  Eigen::VectorXd searchVector(Eigen::Vector3d(0.8, 1, 0));

  std::vector<size_t>             results;
  double                          radius = 0.81; // Two vertices in radius
  bg::model::box<Eigen::VectorXd> search_box(
      searchVector - Eigen::VectorXd::Constant(3, radius),
      searchVector + Eigen::VectorXd::Constant(3, radius));

  tree->query(bg::index::intersects(search_box) and bg::index::satisfies([&](size_t const i) { return bg::distance(searchVector, mesh->vertices().at(i)) <= radius; }),
              std::back_inserter(results));

  BOOST_TEST(results.size() == 2);
  BOOST_TEST(mesh->vertices().at(results.at(0)).getCoords() == Eigen::Vector3d(0, 1, 0));
  BOOST_TEST(mesh->vertices().at(results.at(1)).getCoords() == Eigen::Vector3d(1, 1, 0));
}

/// Resembles how boost geometry is used inside the PetRBF
BOOST_AUTO_TEST_CASE(QueryWithBoxEverything)
{
  PRECICE_TEST(1_rank);
  auto mesh = vertexMesh3D();

  auto tree = rtree::getVertexRTree(mesh);
  BOOST_TEST(tree->size() == 8);

  Eigen::VectorXd searchVector(Eigen::Vector3d(0.8, 1, 0));

  std::vector<size_t>             results;
  double                          radius = std::numeric_limits<double>::max();
  bg::model::box<Eigen::VectorXd> search_box(
      searchVector - Eigen::VectorXd::Constant(3, radius),
      searchVector + Eigen::VectorXd::Constant(3, radius));

  tree->query(bg::index::intersects(search_box) and bg::index::satisfies([&](size_t const i) { return bg::distance(searchVector, mesh->vertices().at(i)) <= radius; }),
              std::back_inserter(results));

  BOOST_TEST(results.size() == 8);
}

BOOST_AUTO_TEST_CASE(QueryWithBoundingBox)
{
  PRECICE_TEST(1_rank);
  auto mesh = vertexMesh3D();

  auto tree = rtree::getVertexRTree(mesh);
  BOOST_TEST(tree->size() == 8);

  Eigen::VectorXd searchVector(Eigen::Vector3d(0.8, 1, 0));

  std::vector<size_t> results;

  auto search_box = toRTreeBox(mesh->getBoundingBox());
  tree->query(bg::index::intersects(search_box), std::back_inserter(results));

  BOOST_TEST(results.size() == tree->size());
}

BOOST_AUTO_TEST_SUITE_END() // Vertex

BOOST_AUTO_TEST_SUITE(Edge)

BOOST_AUTO_TEST_CASE(Query2DEdge)
{
  PRECICE_TEST(1_rank);
  auto mesh = edgeMesh2D();
  auto tree = rtree::getEdgeRTree(mesh);

  BOOST_TEST(tree->size() == 1);

  Eigen::VectorXd     searchVector(Eigen::Vector2d(0.2, 0.8));
  std::vector<size_t> results;

  tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

  BOOST_TEST(results.size() == 1);
  auto &edge = mesh->edges().at(results.front());

  BOOST_TEST(edge.vertex(0).getCoords() == Eigen::Vector2d(1, 0));
  BOOST_TEST(edge.vertex(1).getCoords() == Eigen::Vector2d(1, 1));
}

BOOST_AUTO_TEST_CASE(Query3DEdge)
{
  PRECICE_TEST(1_rank);
  auto mesh = edgeMesh3D();
  auto tree = rtree::getEdgeRTree(mesh);

  BOOST_TEST(tree->size() == 1);

  Eigen::VectorXd     searchVector(Eigen::Vector3d(1.8, 0.0, 0.8));
  std::vector<size_t> results;

  tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

  BOOST_TEST(results.size() == 1);
  auto match = results.front();

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
  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, false, precice::testing::nextMeshID()));
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

  auto tree = rtree::getEdgeRTree(mesh);

  BOOST_TEST(tree->size() == 9);

  Eigen::VectorXd  searchVector(Eigen::Vector3d(0.8, 0.5, 0.0));
  std::set<size_t> results;

  tree->query(bgi::nearest(searchVector, 2), std::inserter(results, results.begin()));

  BOOST_TEST_INFO(results);
  BOOST_TEST(results.size() == 2);
  BOOST_TEST(results.count(eld.getID()) == 1);
  BOOST_TEST(results.count(elr.getID()) == 1);
}

BOOST_AUTO_TEST_SUITE_END() // Edge

BOOST_AUTO_TEST_SUITE(Triangle)

BOOST_AUTO_TEST_CASE(Query3DFullTriangle)
{
  PRECICE_TEST(1_rank);

  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, false, precice::testing::nextMeshID()));
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

  auto tree = rtree::getTriangleRTree(mesh);

  BOOST_TEST(tree->size() == 4);

  Eigen::VectorXd                        searchVector(Eigen::Vector3d(0.7, 0.5, 0.0));
  std::vector<std::pair<double, size_t>> results;

  tree->query(bgi::nearest(searchVector, 3), boost::make_function_output_iterator([&](const precice::mesh::rtree::triangle_traits::IndexType &val) {
                results.push_back(std::make_pair(
                    boost::geometry::distance(
                        searchVector,
                        mesh->triangles().at(val.second)),
                    val.second));
              }));

  std::sort(results.begin(), results.end());
  BOOST_TEST_INFO(results);
  BOOST_TEST(results.size() == 3);
  BOOST_TEST(results.at(0).second == tlb.getID());
  BOOST_TEST(results.at(1).second == tlt.getID());
  BOOST_TEST(results.at(2).second == trt.getID());
  BOOST_TEST(results.at(2).second != trb.getID());
}

BOOST_AUTO_TEST_SUITE_END() // Triangle

BOOST_AUTO_TEST_SUITE(Cache)

BOOST_FIXTURE_TEST_CASE(ClearOnChange, precice::testing::accessors::rtree)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, false, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector2d(0, 0));

  // The Cache should clear whenever a mesh changes
  auto vTree = mesh::rtree::getVertexRTree(mesh);
  BOOST_TEST(getCache().size() == 1);
  mesh->meshChanged(*mesh); // Emit signal, that mesh has changed
  BOOST_TEST(getCache().empty());
}

BOOST_FIXTURE_TEST_CASE(ClearOnDestruction, precice::testing::accessors::rtree)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, false, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector2d(0, 0));

  // The Cache should clear whenever we destroy the Mesh
  auto vTree = mesh::rtree::getVertexRTree(mesh);
  BOOST_TEST(getCache().size() == 1);
  mesh.reset(); // Destroy mesh object, signal is emitted to clear cache
  BOOST_TEST(getCache().empty());
}

BOOST_AUTO_TEST_CASE(CacheVertices)
{
  PRECICE_TEST(1_rank);
  auto ptr = fullMesh();

  auto vt1 = rtree::getVertexRTree(ptr);
  auto vt2 = rtree::getVertexRTree(ptr);
  BOOST_TEST(vt1 == vt2);
}

BOOST_AUTO_TEST_CASE(CacheEdges)
{
  PRECICE_TEST(1_rank);
  auto ptr = fullMesh();

  auto et1 = rtree::getEdgeRTree(ptr);
  auto et2 = rtree::getEdgeRTree(ptr);
  BOOST_TEST(et1 == et2);
}

BOOST_AUTO_TEST_CASE(CacheTriangles)
{
  PRECICE_TEST(1_rank);
  auto ptr = fullMesh();

  auto tt1 = rtree::getTriangleRTree(ptr);
  auto tt2 = rtree::getTriangleRTree(ptr);
  BOOST_TEST(tt1 == tt2);
}

BOOST_AUTO_TEST_CASE(CacheAll)
{
  PRECICE_TEST(1_rank);
  auto ptr = fullMesh();

  auto vt1 = rtree::getVertexRTree(ptr);
  auto et1 = rtree::getEdgeRTree(ptr);
  auto tt1 = rtree::getTriangleRTree(ptr);

  auto vt2 = rtree::getVertexRTree(ptr);
  auto et2 = rtree::getEdgeRTree(ptr);
  auto tt2 = rtree::getTriangleRTree(ptr);

  BOOST_TEST(vt1 == vt2);
  BOOST_TEST(et1 == et2);
  BOOST_TEST(tt1 == tt2);
}

BOOST_AUTO_TEST_SUITE_END() // Cache

BOOST_AUTO_TEST_SUITE_END() // Query
BOOST_AUTO_TEST_SUITE_END() // RTree
BOOST_AUTO_TEST_SUITE_END() // Mesh
