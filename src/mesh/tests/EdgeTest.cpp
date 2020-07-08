#include <Eigen/Core>
#include <iosfwd>
#include <string>
#include "logging/Logger.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(EdgeTests)

BOOST_AUTO_TEST_CASE(Edges)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Eigen::Vector3d::Constant(0.0), 0);
  Vertex v2(Eigen::Vector3d::Constant(1.0), 1);
  Edge   edge(v1, v2, 0);

  Eigen::VectorXd coords1 = edge.vertex(0).getCoords();
  Eigen::VectorXd coords2 = edge.vertex(1).getCoords();
  BOOST_TEST(coords1 == Eigen::Vector3d::Constant(0.0));
  BOOST_TEST(coords2 == Eigen::Vector3d::Constant(1.0));
}

BOOST_AUTO_TEST_CASE(EdgeEquality)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Eigen::Vector3d(0, 0, 0), 0);
  Vertex v2(Eigen::Vector3d(0, 0, 1), 0);
  Vertex v3(Eigen::Vector3d(0, 0, 2), 0);
  Edge   edge1(v1, v2, 0);
  Edge   edge2(v2, v1, 1);
  Edge   edge3(v1, v3, 0);
  Edge   edge4(v1, v3, 0);
  edge4.setNormal(Eigen::Vector3d(Eigen::Vector3d(0, 1, 0)));
  BOOST_TEST(edge1 == edge2);
  BOOST_TEST(edge1 != edge3);
  BOOST_TEST(edge3 != edge4);
}

BOOST_AUTO_TEST_CASE(EdgeWKTPrint)
{
  PRECICE_TEST(1_rank);
  Vertex            v1(Eigen::Vector2d(1., 2.), 0);
  Vertex            v2(Eigen::Vector2d(2., 3.), 0);
  Edge              e1(v1, v2, 0);
  std::stringstream e1stream;
  e1stream << e1;
  std::string e1str("LINESTRING (1 2, 2 3)");
  BOOST_TEST(e1str == e1stream.str());
  Vertex            v3(Eigen::Vector3d(1., 2., 3.), 0);
  Vertex            v4(Eigen::Vector3d(3., 2., 1.), 0);
  Edge              e2(v3, v4, 0);
  std::stringstream e2stream;
  e2stream << e2;
  std::string e2str("LINESTRING (1 2 3, 3 2 1)");
  BOOST_TEST(e2str == e2stream.str());
}

BOOST_AUTO_TEST_CASE(EdgeConnectedTo)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Eigen::Vector3d(0, 0, 1), 0);
  Vertex v2(Eigen::Vector3d(0, 0, 2), 0);
  Vertex v3(Eigen::Vector3d(0, 0, 3), 0);
  Vertex v4(Eigen::Vector3d(0, 0, 4), 0);

  Edge edge1(v1, v2, 0);
  Edge edge2(v2, v3, 0);
  BOOST_TEST(edge1.connectedTo(edge2));
  BOOST_TEST(edge2.connectedTo(edge1));

  Edge edge3(v3, v4, 0);
  BOOST_TEST(!edge1.connectedTo(edge3));
  BOOST_TEST(!edge3.connectedTo(edge1));
  BOOST_TEST(edge2.connectedTo(edge3));
  BOOST_TEST(edge3.connectedTo(edge2));
}

BOOST_AUTO_TEST_SUITE_END() // Edge
BOOST_AUTO_TEST_SUITE_END() // Mesh
