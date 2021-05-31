#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <iosfwd>
#include <string>
#include "logging/Logger.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(EdgeTests)

BOOST_AUTO_TEST_CASE(Edges)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector3d::Constant(0.0), 0);
  Vertex v2(Vector3d::Constant(1.0), 1);
  Edge   edge(v1, v2, 0);

  BOOST_TEST(edge.getDimensions() == 3);
  VectorXd coords1 = edge.vertex(0).getCoords();
  VectorXd coords2 = edge.vertex(1).getCoords();
  BOOST_TEST(coords1 == Vector3d::Constant(0.0));
  BOOST_TEST(coords2 == Vector3d::Constant(1.0));

  double expectedLenght = std::sqrt(3.0);
  double expectedRadius = expectedLenght / 2.0;
  BOOST_TEST(edge.getLength() == expectedLenght);
  BOOST_TEST(edge.getEnclosingRadius() == expectedRadius);
  BOOST_TEST(testing::equals(edge.getCenter(), (coords1 + coords2) / 2));
  BOOST_TEST(edge.computeNormal().dot(coords1 - coords2) == 0.0);
}

BOOST_AUTO_TEST_CASE(Dimensions2D)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector2d::Constant(0.0), 0);
  Vertex v2(Vector2d::Constant(1.0), 1);
  Edge   edge(v1, v2, 0);
  BOOST_TEST(edge.getDimensions() == 2);

  double expectedLenght = std::sqrt(2.0);
  double expectedRadius = expectedLenght / 2.0;
  BOOST_TEST(edge.getLength() == expectedLenght);
  BOOST_TEST(edge.getEnclosingRadius() == expectedRadius);
  BOOST_TEST(testing::equals(edge.getCenter(), Vector2d{0.5, 0.5}));
  BOOST_TEST(edge.computeNormal().dot(Vector2d{0.5, 0.5}) == 0.0);
}

BOOST_AUTO_TEST_CASE(Dimensions2DX)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector2d::Constant(0.0), 0);
  Vertex v2(Vector2d{1, 0}, 1);
  Edge   edge(v1, v2, 0);
  BOOST_TEST(edge.getDimensions() == 2);

  double expectedLenght = std::sqrt(1.0);
  double expectedRadius = expectedLenght / 2.0;
  BOOST_TEST(edge.getLength() == expectedLenght);
  BOOST_TEST(edge.getEnclosingRadius() == expectedRadius);
  BOOST_TEST(testing::equals(edge.getCenter(), Vector2d{0.5, 0}));
  BOOST_TEST(edge.computeNormal().dot(Vector2d{1, 0}) == 0.0);
}

BOOST_AUTO_TEST_CASE(Dimensions2DY)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector2d::Constant(0.0), 0);
  Vertex v2(Vector2d{0, 1}, 1);
  Edge   edge(v1, v2, 0);
  BOOST_TEST(edge.getDimensions() == 2);

  double expectedLenght = std::sqrt(1.0);
  double expectedRadius = expectedLenght / 2.0;
  BOOST_TEST(edge.getLength() == expectedLenght);
  BOOST_TEST(edge.getEnclosingRadius() == expectedRadius);
  BOOST_TEST(testing::equals(edge.getCenter(), Vector2d{0, 0.5}));
  BOOST_TEST(edge.computeNormal().dot(Vector2d{0, 1}) == 0.0);
}

BOOST_AUTO_TEST_CASE(Dimensions3DX)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector3d::Constant(0.0), 0);
  Vertex v2(Vector3d{1, 0, 0}, 1);
  Edge   edge(v1, v2, 0);
  BOOST_TEST(edge.getDimensions() == 3);

  double expectedLenght = std::sqrt(1.0);
  double expectedRadius = expectedLenght / 2.0;
  BOOST_TEST(edge.getLength() == expectedLenght);
  BOOST_TEST(edge.getEnclosingRadius() == expectedRadius);
  BOOST_TEST(testing::equals(edge.getCenter(), Vector3d{0.5, 0, 0}));
  const Vector3d normal = edge.computeNormal();
  BOOST_TEST(normal.norm() == 1.0);
  BOOST_TEST(normal.dot(Vector3d{1, 0, 0}) == 0.0);
}

BOOST_AUTO_TEST_CASE(Dimensions3DY)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector3d::Constant(0.0), 0);
  Vertex v2(Vector3d{0, 1, 0}, 1);
  Edge   edge(v1, v2, 0);
  BOOST_TEST(edge.getDimensions() == 3);

  double expectedLenght = std::sqrt(1.0);
  double expectedRadius = expectedLenght / 2.0;
  BOOST_TEST(edge.getLength() == expectedLenght);
  BOOST_TEST(edge.getEnclosingRadius() == expectedRadius);
  BOOST_TEST(testing::equals(edge.getCenter(), Vector3d{0, 0.5, 0}));
  const Vector3d normal = edge.computeNormal();
  BOOST_TEST(normal.norm() == 1.0);
  BOOST_TEST(normal.dot(Vector3d{0, 1, 0}) == 0.0);
}

BOOST_AUTO_TEST_CASE(Dimensions3DZ)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector3d::Constant(0.0), 0);
  Vertex v2(Vector3d{0, 0, 1}, 1);
  Edge   edge(v1, v2, 0);
  BOOST_TEST(edge.getDimensions() == 3);

  double expectedLenght = std::sqrt(1.0);
  double expectedRadius = expectedLenght / 2.0;
  BOOST_TEST(edge.getLength() == expectedLenght);
  BOOST_TEST(edge.getEnclosingRadius() == expectedRadius);
  BOOST_TEST(testing::equals(edge.getCenter(), Vector3d{0, 0, 0.5}));
  const Vector3d normal = edge.computeNormal();
  BOOST_TEST(normal.norm() == 1.0);
  BOOST_TEST(normal.dot(Vector3d{0, 0, 1}) == 0.0);
}

BOOST_AUTO_TEST_CASE(Dimensions3D)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector3d::Constant(0.0), 0);
  Vertex v2(Vector3d::Constant(1.0), 1);
  Edge   edge(v1, v2, 0);
  BOOST_TEST(edge.getDimensions() == 3);

  double expectedLenght = std::sqrt(3.0);
  double expectedRadius = expectedLenght / 2.0;
  BOOST_TEST(edge.getLength() == expectedLenght);
  BOOST_TEST(edge.getEnclosingRadius() == expectedRadius);
  BOOST_TEST(testing::equals(edge.getCenter(), Vector3d::Constant(0.5)));
  BOOST_TEST(edge.computeNormal().dot(Vector3d::Constant(0.3333)) == 0.0);
}

BOOST_AUTO_TEST_CASE(EdgeEquality)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector3d(0, 0, 0), 0);
  Vertex v2(Vector3d(0, 0, 1), 0);
  Vertex v3(Vector3d(0, 0, 2), 0);
  Edge   edge1(v1, v2, 0);
  Edge   edge2(v2, v1, 1);
  Edge   edge3(v1, v3, 0);
  Edge   edge4(v1, v3, 0);
  BOOST_TEST(edge1 == edge2);
  BOOST_TEST(edge1 != edge3);
  BOOST_TEST(edge3 == edge4);
}

BOOST_AUTO_TEST_CASE(ComputeNormal2D_Unit)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector2d{0.0, 0.0}, 0);
  Vertex v2(Vector2d{1.0, 0.0}, 1);
  Edge   edge(v1, v2, 0);

  auto normal = edge.computeNormal();
  BOOST_TEST(normal.size() == 2);
  BOOST_TEST(normal.norm() == 1.0);
  Vector2d a{0.0, 1.0};
  Vector2d b{0.0, -1.0};
  BOOST_TEST((normal == a || normal == b));
}

BOOST_AUTO_TEST_CASE(ComputeNormal2D)
{
  PRECICE_TEST(1_rank);
  Vector2d a{-0.5, 1.0};
  Vector2d b{1.25, -1.1};
  Vertex   v1(a, 0);
  Vertex   v2(b, 1);
  Edge     edge(v1, v2, 0);

  auto normal = edge.computeNormal();
  BOOST_TEST(normal.size() == 2);
  BOOST_TEST(normal.norm() == 1.0);
  BOOST_TEST(normal.dot(b - a) == 0.0);
}

BOOST_AUTO_TEST_CASE(ComputeNormal3D_Unit)
{
  PRECICE_TEST(1_rank);
  Vector3d a{0.0, 0.0, 1.0};
  Vector3d b{1.0, 0.0, 1.0};
  Vertex   v1(a, 0);
  Vertex   v2(b, 1);
  Edge     edge(v1, v2, 0);

  auto normal = edge.computeNormal();
  BOOST_TEST(normal.size() == 3);
  BOOST_TEST(normal.norm() == 1.0);
  BOOST_TEST(normal.dot(b - a) == 0.0);
}

BOOST_AUTO_TEST_CASE(EdgeWKTPrint)
{
  PRECICE_TEST(1_rank);
  Vertex            v1(Vector2d(1., 2.), 0);
  Vertex            v2(Vector2d(2., 3.), 0);
  Edge              e1(v1, v2, 0);
  std::stringstream e1stream;
  e1stream << e1;
  std::string e1str("LINESTRING (1 2, 2 3)");
  BOOST_TEST(e1str == e1stream.str());
  Vertex            v3(Vector3d(1., 2., 3.), 0);
  Vertex            v4(Vector3d(3., 2., 1.), 0);
  Edge              e2(v3, v4, 0);
  std::stringstream e2stream;
  e2stream << e2;
  std::string e2str("LINESTRING (1 2 3, 3 2 1)");
  BOOST_TEST(e2str == e2stream.str());
}

BOOST_AUTO_TEST_CASE(EdgeConnectedTo)
{
  PRECICE_TEST(1_rank);
  Vertex v1(Vector3d(0, 0, 1), 0);
  Vertex v2(Vector3d(0, 0, 2), 0);
  Vertex v3(Vector3d(0, 0, 3), 0);
  Vertex v4(Vector3d(0, 0, 4), 0);

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
