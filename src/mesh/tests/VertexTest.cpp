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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Vertices)
{
  PRECICE_TEST();
  mesh::Vertex vertex(Eigen::Vector3d::Constant(1.0), 0);

  Eigen::Vector3d coords = vertex.getCoords();
  BOOST_TEST(testing::equals(coords, Eigen::Vector3d::Constant(1.0)));

  int id = vertex.getID();
  BOOST_TEST(id == 0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VertexEquality)
{
  PRECICE_TEST();
  using namespace mesh;
  using namespace Eigen;
  Vertex v1(Vector3d::Constant(4.0), 0);
  Vertex v2(Vector3d::Constant(4.0), 1);
  Vertex v3(Vector3d::Constant(2.0), 0);
  BOOST_TEST(v1 == v2);
  BOOST_TEST(v1 != v3);
  BOOST_TEST(v2 != v3);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VertexWKTPrint)
{
  PRECICE_TEST();
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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(RawCoordsNorm)
{
  PRECICE_TEST();
  using namespace mesh;
  // 2D vertices are padded with zeros, so the third component cancels
  Vertex a2(Eigen::Vector2d(1., 2.), 0);
  Vertex b2(Eigen::Vector2d(4., 6.), 1);
  BOOST_TEST(norm(a2.rawCoords(), b2.rawCoords()) == 5.0);
  BOOST_TEST(norm(a2.rawCoords(), b2.rawCoords()) == (a2.getCoords() - b2.getCoords()).norm());

  Vertex a3(Eigen::Vector3d(1., -2., 3.), 0);
  Vertex b3(Eigen::Vector3d(-0.5, 4., 0.25), 1);
  BOOST_TEST(norm(a3.rawCoords(), b3.rawCoords()) == (a3.getCoords() - b3.getCoords()).norm(), boost::test_tools::tolerance(1e-14));
  BOOST_TEST(norm(a3.rawCoords(), b3.rawCoords()) == norm(b3.rawCoords(), a3.rawCoords()));
  BOOST_TEST(norm(a3.rawCoords(), a3.rawCoords()) == 0.0);
}

BOOST_AUTO_TEST_SUITE_END() // Vertex
BOOST_AUTO_TEST_SUITE_END() // Mesh
