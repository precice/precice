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
}

BOOST_AUTO_TEST_CASE(VertexEquality)
{
    using namespace mesh;
    using namespace Eigen;
    Vertex v1(Vector3d::Constant(4.0), 0);
    Vertex v2(Vector3d::Constant(4.0), 1);
    Vertex v3(Vector3d::Constant(2.0), 0);
    BOOST_TEST(v1 == v2);
    BOOST_TEST(v1 != v3);
    BOOST_TEST(v2 != v3);
}
BOOST_AUTO_TEST_SUITE_END() // Mesh
