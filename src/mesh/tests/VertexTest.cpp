#include <Eigen/Core>
#include <iosfwd>
#include <string>
#include "logging/Logger.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(VertexTests)

BOOST_AUTO_TEST_CASE(Vertices)
{
  PRECICE_TEST(1_rank);
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
  PRECICE_TEST(1_rank);
  using namespace mesh;
  using namespace Eigen;
  Vertex v1(Vector3d::Constant(4.0), 0);
  Vertex v2(Vector3d::Constant(4.0), 1);
  Vertex v3(Vector3d::Constant(2.0), 0);
  BOOST_TEST(v1 == v2);
  BOOST_TEST(v1 != v3);
  BOOST_TEST(v2 != v3);
}

BOOST_AUTO_TEST_CASE(VertexWKTPrint)
{
  PRECICE_TEST(1_rank);
  using namespace mesh;
  Vertex            v1(Eigen::Vector2d(1., 2.), 0);
  std::stringstream v1stream;
  v1stream << v1;
  std::string v1str("POINT (1 2)");
  BOOST_TEST(v1str == v1stream.str());
  Vertex            v2(Eigen::Vector3d(1., 2., 3.), 0);
  std::stringstream v2stream;
  v2stream << v2;
  std::string v2str("POINT (1 2 3)");
  BOOST_TEST(v2str == v2stream.str());
}

BOOST_AUTO_TEST_SUITE_END() // Vertex
BOOST_AUTO_TEST_SUITE_END() // Mesh
