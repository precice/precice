#include "testing/Testing.hpp"
#include "mesh/RTree.hpp"
#include "mesh/impl/RTreeAdapter.hpp"
#include "math/geometry.hpp"

using namespace precice::mesh;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(RTree)

BOOST_AUTO_TEST_CASE(Query_2D)
{
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, false));
  mesh->createVertex(Eigen::Vector2d(0, 0));
  mesh->createVertex(Eigen::Vector2d(0, 1));
  mesh->createVertex(Eigen::Vector2d(1, 0));
  mesh->createVertex(Eigen::Vector2d(1, 1));
  
  auto tree = rtree::getVertexRTree(mesh);

  BOOST_TEST(tree->size() == 4);

  Eigen::VectorXd searchVector(Eigen::Vector2d(0.2, 0.8));
  std::vector<size_t> results;

  tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

  BOOST_TEST(results.size() == 1);
  BOOST_TEST( mesh->vertices()[results[0]].getCoords() == Eigen::Vector2d(0, 1) );
}

BOOST_AUTO_TEST_CASE(Query_3D)
{
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 3, false));
  mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  mesh->createVertex(Eigen::Vector3d(0, 0, 1));
  mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  mesh->createVertex(Eigen::Vector3d(0, 1, 1));
  mesh->createVertex(Eigen::Vector3d(1, 0, 0));
  mesh->createVertex(Eigen::Vector3d(1, 0, 1));
  mesh->createVertex(Eigen::Vector3d(1, 1, 0));
  mesh->createVertex(Eigen::Vector3d(1, 1, 1));

  auto tree = rtree::getVertexRTree(mesh);

  BOOST_TEST(tree->size() == 8);

  Eigen::VectorXd searchVector(Eigen::Vector3d(0.8, 0.0, 0.8));
  std::vector<size_t> results;

  tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

  BOOST_TEST(results.size() == 1);
  BOOST_TEST( mesh->vertices()[results[0]].getCoords() == Eigen::Vector3d(1, 0, 1) );
}

/// Resembles how boost geometry is used inside the PetRBF
BOOST_AUTO_TEST_CASE(QueryWithBox)
{
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 3, false));
  mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  mesh->createVertex(Eigen::Vector3d(0, 0, 1));
  mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  mesh->createVertex(Eigen::Vector3d(0, 1, 1));
  mesh->createVertex(Eigen::Vector3d(1, 0, 0));
  mesh->createVertex(Eigen::Vector3d(1, 0, 1));
  mesh->createVertex(Eigen::Vector3d(1, 1, 0));
  mesh->createVertex(Eigen::Vector3d(1, 1, 1));

  auto tree = rtree::getVertexRTree(mesh);
  BOOST_TEST(tree->size() == 8);
  
  Eigen::VectorXd searchVector(Eigen::Vector3d(0.8, 1, 0));

  {
    std::vector<size_t> results;
    double radius = 0.1; // No vertices in radius
    bg::model::box<Eigen::VectorXd> search_box(
      searchVector - Eigen::VectorXd::Constant(3, radius),
      searchVector + Eigen::VectorXd::Constant(3, radius));
    
    tree->query(bg::index::within(search_box) and bg::index::satisfies([&](size_t const i){
          return bg::distance(searchVector, mesh->vertices()[i]) <= radius;}),
      std::back_inserter(results));
  
    BOOST_TEST(results.size() == 0);
  }

  {
    std::vector<size_t> results;
    double radius = 0.81; // Two vertices in radius
    bg::model::box<Eigen::VectorXd> search_box(
      searchVector - Eigen::VectorXd::Constant(3, radius),
      searchVector + Eigen::VectorXd::Constant(3, radius));
    
    tree->query(bg::index::within(search_box) and bg::index::satisfies([&](size_t const i){
          return bg::distance(searchVector, mesh->vertices()[i]) <= radius;}),
      std::back_inserter(results));
    
    BOOST_TEST(results.size() == 2);
    BOOST_TEST(mesh->vertices()[results[0]].getCoords() == Eigen::Vector3d(0, 1, 0));
    BOOST_TEST(mesh->vertices()[results[1]].getCoords() == Eigen::Vector3d(1, 1, 0));
  }

  {
    std::vector<size_t> results;
    double radius = std::numeric_limits<double>::max();
    bg::model::box<Eigen::VectorXd> search_box(
      searchVector - Eigen::VectorXd::Constant(3, radius),
      searchVector + Eigen::VectorXd::Constant(3, radius));
    
    tree->query(bg::index::within(search_box) and bg::index::satisfies([&](size_t const i){
          return bg::distance(searchVector, mesh->vertices()[i]) <= radius;}),
      std::back_inserter(results));
    
    BOOST_TEST(results.size() == 8);
  }
}

BOOST_AUTO_TEST_CASE(VertexAdapter)
{
  precice::mesh::Mesh mesh("MyMesh", 2, false);
  auto & v = mesh.createVertex(Eigen::Vector2d(1, 2));
  BOOST_TEST(bg::get<0>(v) == 1);
  BOOST_TEST(bg::get<1>(v) == 2);
  BOOST_TEST(bg::get<2>(v) == 0);
  bg::set<1>(v, 5);
  BOOST_TEST(bg::get<1>(v) == 5);
}

BOOST_AUTO_TEST_CASE(EdgeAdapter)
{
  precice::mesh::Mesh mesh("MyMesh", 2, false);
  auto & v1 = mesh.createVertex(Eigen::Vector2d(1, 2));
  auto & v2 = mesh.createVertex(Eigen::Vector2d(3, 4));
  auto & e = mesh.createEdge(v1, v2);
  BOOST_TEST((bg::get<0,0>(e)) == 1.0);
  BOOST_TEST((bg::get<0,1>(e)) == 2.0);
  BOOST_TEST((bg::get<0,2>(e)) == 0.0);
  BOOST_TEST((bg::get<1,0>(e)) == 3.0);
  BOOST_TEST((bg::get<1,1>(e)) == 4.0);
  BOOST_TEST((bg::get<1,2>(e)) == 0.0);
  bg::set<1,1>(e, 5.0);
  BOOST_TEST((bg::get<1,1>(e)) == 5.0);
}

BOOST_AUTO_TEST_CASE(TriangleAdapter)
{
  precice::mesh::Mesh mesh("MyMesh", 3, false);
  auto & v1 = mesh.createVertex(Eigen::Vector3d(0, 2, 0));
  auto & v2 = mesh.createVertex(Eigen::Vector3d(2, 1, 0));
  auto & v3 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto & e1 = mesh.createEdge(v1, v2);
  auto & e2 = mesh.createEdge(v2, v3);
  auto & e3 = mesh.createEdge(v3, v1);
  auto & t = mesh.createTriangle(e1, e2, e3);

  std::vector<Eigen::VectorXd> vertices(t.begin(), t.end());
  std::vector<Eigen::VectorXd> refs{ v1.getCoords(), v2.getCoords(), v3.getCoords()};
  BOOST_TEST(vertices.size() == refs.size());
  BOOST_TEST((std::is_permutation(
                  vertices.begin(), vertices.end(),
                  refs.begin(),
                  []( const Eigen::VectorXd& lhs, const Eigen::VectorXd& rhs) {
                     return precice::math::equals(lhs, rhs);
                  }
            )));
}

BOOST_AUTO_TEST_CASE(QuadAdapter)
{
  precice::mesh::Mesh mesh("MyMesh", 3, false);
  auto & v1 = mesh.createVertex(Eigen::Vector3d(0, 2, 0));
  auto & v2 = mesh.createVertex(Eigen::Vector3d(2, 1, 0));
  auto & v3 = mesh.createVertex(Eigen::Vector3d(3, 0, 0));
  auto & v4 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto & e1 = mesh.createEdge(v1, v2);
  auto & e2 = mesh.createEdge(v2, v3);
  auto & e3 = mesh.createEdge(v3, v4);
  auto & e4 = mesh.createEdge(v4, v1);
  auto & t = mesh.createQuad(e1, e2, e3, e4);

  std::vector<Eigen::VectorXd> vertices(t.begin(), t.end());
  std::vector<Eigen::VectorXd> refs{ v1.getCoords(), v2.getCoords(), v3.getCoords(), v4.getCoords()};
  BOOST_TEST(vertices.size() == refs.size());
  BOOST_TEST((std::is_permutation(
                  vertices.begin(), vertices.end(),
                  refs.begin(),
                  []( const Eigen::VectorXd& lhs, const Eigen::VectorXd& rhs) {
                     return precice::math::equals(lhs, rhs);
                  }
            )));
}

BOOST_AUTO_TEST_CASE(VectorAdapter)
{
  Eigen::VectorXd vec = Eigen::Vector2d(1, 2);
  BOOST_TEST(bg::get<0>(vec) == 1);
  BOOST_TEST(bg::get<1>(vec) == 2);
  BOOST_TEST(bg::get<2>(vec) == 0);
  bg::set<1>(vec, 5);
  BOOST_TEST(bg::get<1>(vec) == 5);
}

BOOST_AUTO_TEST_CASE(CacheClearing)
{
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, false));
  mesh->createVertex(Eigen::Vector2d(0, 0));
  
  auto tree1 = rtree::getVertexRTree(mesh);
  BOOST_TEST(rtree::trees.size() == 1);
  mesh->meshChanged(*mesh); // Emit signal, that mesh has changed
  BOOST_TEST(rtree::trees.size() == 0);
  
  auto tree2 = rtree::getVertexRTree(mesh);
  BOOST_TEST(rtree::trees.size() == 1);
  mesh.reset(); // Destroy mesh object, signal is emitted to clear cache
  BOOST_TEST(rtree::trees.size() == 0);
  
}


BOOST_AUTO_TEST_SUITE_END() // RTree
BOOST_AUTO_TEST_SUITE_END() // Mesh
