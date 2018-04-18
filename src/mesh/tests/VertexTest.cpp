#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"
#include <Eigen/Core>

using namespace precice;

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(Vertices)
{
  mesh::Vertex vertex(Eigen::Vector3d::Constant(1.0), 0);

  Eigen::Vector3d coords = vertex.getCoords();
  BOOST_TEST(testing::equals(coords, Eigen::Vector3d::Constant(1.0)));

  int id = vertex.getID();
  BOOST_TEST(id == 0);

  Eigen::Vector3d normal = vertex.getNormal();
  BOOST_TEST(testing::equals(normal, Eigen::Vector3d::Zero()));

  void *mesh = static_cast<void *>(vertex.mesh());
  BOOST_TEST(mesh == static_cast<void *>(nullptr));
}

BOOST_AUTO_TEST_SUITE_END() // Mesh
