#include "testing/Testing.hpp"
#include "mesh/RTree.hpp"
#include "mesh/impl/RTreeAdapter.hpp"

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

BOOST_AUTO_TEST_CASE(VectorAdapter)
{
  Eigen::VectorXd vec = Eigen::Vector2d(1, 2);
  BOOST_TEST(bg::get<0>(vec) == 1);
  BOOST_TEST(bg::get<1>(vec) == 2);
  BOOST_TEST(bg::get<2>(vec) == 0);
  bg::set<1>(vec, 5);
  BOOST_TEST(bg::get<1>(vec) == 5);
}


BOOST_AUTO_TEST_SUITE_END() // RTree
BOOST_AUTO_TEST_SUITE_END() // Mesh
