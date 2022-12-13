#include <vector>
#include "Eigen/Core"
#include "mapping/Polation.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mapping;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(Interpolation)

BOOST_AUTO_TEST_CASE(VertexInterpolation)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d location(0.0, 0.0, 0.0);
  mesh::Vertex    vertex(Eigen::Vector3d(1.0, 2.0, 0.0), 0);

  Polation polation(location, vertex);

  std::vector<int>    expectedIndices = {0};
  std::vector<double> expectedWeights = {1.0};

  BOOST_TEST(polation.getWeightedElements().size() == 1);
  BOOST_TEST(polation.isInterpolation());
  BOOST_TEST(polation.distance() == vertex.getCoords().norm());

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {
    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_CASE(EdgeInterpolation)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex v1(Eigen::Vector3d(0.0, 0.0, 0.0), 0);
  mesh::Vertex v2(Eigen::Vector3d(0.0, 2.0, 0.0), 1);
  mesh::Edge   edge(v1, v2);

  Eigen::Vector3d location(0.0, 0.4, 0.0);

  Polation polation(location, edge);

  std::vector<int>    expectedIndices = {0, 1};
  std::vector<double> expectedWeights = {0.8, 0.2};

  BOOST_TEST(polation.getWeightedElements().size() == 2);
  BOOST_TEST(polation.isInterpolation());
  BOOST_TEST(polation.distance() == 0);

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {
    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_CASE(EdgeProjectedInterpolation)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex v1(Eigen::Vector3d(0.0, 0.0, 0.0), 0);
  mesh::Vertex v2(Eigen::Vector3d(0.0, 2.0, 0.0), 1);
  mesh::Edge   edge(v1, v2);

  Eigen::Vector3d location(0.0, 0.4, 0.12);

  Polation polation(location, edge);

  std::vector<int>    expectedIndices = {0, 1};
  std::vector<double> expectedWeights = {0.8, 0.2};

  BOOST_TEST(polation.getWeightedElements().size() == 2);
  BOOST_TEST(polation.isInterpolation());
  BOOST_TEST(polation.distance() == 0.12);

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {
    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_CASE(TriangleInterpolation)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex   v1(Eigen::Vector3d(0.0, 0.0, 0.0), 0);
  mesh::Vertex   v2(Eigen::Vector3d(2.0, 0.0, 0.0), 1);
  mesh::Vertex   v3(Eigen::Vector3d(1.0, 2.0, 0.0), 2);
  mesh::Edge     e1(v1, v2);
  mesh::Edge     e2(v2, v3);
  mesh::Edge     e3(v1, v3);
  mesh::Triangle triangle(e1, e2, e3);
  triangle.computeNormal();

  Eigen::Vector3d location(1.0, 0.6, 0.0);

  Polation polation(location, triangle);

  std::vector<int>    expectedIndices = {0, 1, 2};
  std::vector<double> expectedWeights = {0.35, 0.35, 0.3};

  BOOST_TEST(polation.getWeightedElements().size() == 3);
  BOOST_TEST(polation.isInterpolation());
  BOOST_TEST(polation.distance() == 0);

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {

    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_CASE(TriangleProjectedInterpolation)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex   v1(Eigen::Vector3d(0.0, 0.0, 0.0), 0);
  mesh::Vertex   v2(Eigen::Vector3d(2.0, 0.0, 0.0), 1);
  mesh::Vertex   v3(Eigen::Vector3d(1.0, 2.0, 0.0), 2);
  mesh::Edge     e1(v1, v2);
  mesh::Edge     e2(v2, v3);
  mesh::Edge     e3(v1, v3);
  mesh::Triangle triangle(e1, e2, e3);
  triangle.computeNormal();

  Eigen::Vector3d location(1.0, 0.6, 0.14);

  Polation polation(location, triangle);

  std::vector<int>    expectedIndices = {0, 1, 2};
  std::vector<double> expectedWeights = {0.35, 0.35, 0.3};

  BOOST_TEST(polation.getWeightedElements().size() == 3);
  BOOST_TEST(polation.isInterpolation());
  BOOST_TEST(polation.distance() == 0.14);

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {

    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_CASE(EdgeExtrapolation)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex v1(Eigen::Vector3d(0.0, 0.0, 0.0), 0);
  mesh::Vertex v2(Eigen::Vector3d(0.0, 2.0, 0.0), 1);
  mesh::Edge   edge(v1, v2);

  Eigen::Vector3d location(0.0, 3.0, 0.0);

  Polation polation(location, edge);

  std::vector<int>    expectedIndices = {0, 1};
  std::vector<double> expectedWeights = {-0.5, 1.5};

  BOOST_TEST(polation.getWeightedElements().size() == 2);
  BOOST_TEST(not polation.isInterpolation());
  // The distance to projection is still 0
  BOOST_TEST(polation.distance() == 0);

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {
    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_CASE(TriangleExtrapolation)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex   v1(Eigen::Vector3d(0.0, 0.0, 0.0), 0);
  mesh::Vertex   v2(Eigen::Vector3d(2.0, 0.0, 0.0), 1);
  mesh::Vertex   v3(Eigen::Vector3d(1.0, 2.0, 0.0), 2);
  mesh::Edge     e1(v1, v2);
  mesh::Edge     e2(v2, v3);
  mesh::Edge     e3(v1, v3);
  mesh::Triangle triangle(e1, e2, e3);
  triangle.computeNormal();

  Eigen::Vector3d location(4.0, 0.6, 0.0);

  Polation polation(location, triangle);

  std::vector<int>    expectedIndices = {0, 1, 2};
  std::vector<double> expectedWeights = {-1.15, 1.85, 0.3};

  BOOST_TEST(polation.getWeightedElements().size() == 3);
  BOOST_TEST(not polation.isInterpolation());
  // The distance to projection is still 0
  BOOST_TEST(polation.distance() == 0);

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {

    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_CASE(TetrahedronInterpolation)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex v1(Eigen::Vector3d(1.0, 0.0, 0.0), 0);
  mesh::Vertex v2(Eigen::Vector3d(0.0, 1.0, 0.0), 1);
  mesh::Vertex v3(Eigen::Vector3d(0.0, 0.0, 1.0), 2);
  mesh::Vertex v4(Eigen::Vector3d(0.0, 0.0, 0.0), 3);

  mesh::Tetrahedron tetra(v1, v2, v3, v4);

  Eigen::Vector3d location(0.15, 0.25, 0.40);

  Polation polation(location, tetra);

  std::vector<int>    expectedIndices = {0, 1, 2, 3};
  std::vector<double> expectedWeights = {0.15, 0.25, 0.40, 0.20};

  BOOST_TEST(polation.getWeightedElements().size() == 4);
  BOOST_TEST(polation.isInterpolation());
  BOOST_TEST(polation.distance() == 0);

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {

    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_CASE(TetrahedronExtrapolation)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex v1(Eigen::Vector3d(1.0, 0.0, 0.0), 0);
  mesh::Vertex v2(Eigen::Vector3d(0.0, 1.0, 0.0), 1);
  mesh::Vertex v3(Eigen::Vector3d(0.0, 0.0, 1.0), 2);
  mesh::Vertex v4(Eigen::Vector3d(0.0, 0.0, 0.0), 3);

  mesh::Tetrahedron tetra(v1, v2, v3, v4);

  Eigen::Vector3d location(-0.15, 0.25, 0.40);

  Polation polation(location, tetra);

  std::vector<int>    expectedIndices = {0, 1, 2, 3};
  std::vector<double> expectedWeights = {-0.15, 0.25, 0.40, 0.50};

  BOOST_TEST(polation.getWeightedElements().size() == 4);
  BOOST_TEST(not polation.isInterpolation());
  // There is no projection distance as such. This should always be 0
  BOOST_TEST(polation.distance() == 0);

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {

    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_CASE(PolationToleranceEdge)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex v1(Eigen::Vector3d(1.0, 0.0, 0.0), 0);
  mesh::Vertex v2(Eigen::Vector3d(0.0, 1.0, 0.0), 1);

  mesh::Edge edge(v1, v2);

  Eigen::Vector3d location(1 + 1e-15, -1e-16, 0.0);

  Polation polation(location, edge);

  BOOST_TEST(polation.isInterpolation());
}

BOOST_AUTO_TEST_CASE(PolationToleranceTriangle)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex v1(Eigen::Vector3d(1.0, 0.0, 0.0), 0);
  mesh::Vertex v2(Eigen::Vector3d(0.0, 1.0, 0.0), 1);
  mesh::Vertex v3(Eigen::Vector3d(0.0, 0.0, 0.0), 2);

  mesh::Triangle triangle(v1, v2, v3);

  Eigen::Vector3d location(0.5 + 1e-15, 0.5 + 1e-15, 1e-14);

  Polation polation(location, triangle);

  BOOST_TEST(polation.isInterpolation());
}

BOOST_AUTO_TEST_CASE(PolationToleranceTetra)
{
  PRECICE_TEST(1_rank);
  mesh::Vertex v1(Eigen::Vector3d(1.0, 0.0, 0.0), 0);
  mesh::Vertex v2(Eigen::Vector3d(0.0, 1.0, 0.0), 1);
  mesh::Vertex v3(Eigen::Vector3d(0.0, 0.0, 1.0), 2);
  mesh::Vertex v4(Eigen::Vector3d(0.0, 0.0, 0.0), 3);

  mesh::Tetrahedron tetra(v1, v2, v3, v4);

  Eigen::Vector3d location(1 - 1e-15, -1e-15, -1e-15);

  Polation polation(location, tetra);

  BOOST_TEST(polation.isInterpolation());
}

BOOST_AUTO_TEST_SUITE_END() // Interpolation
BOOST_AUTO_TEST_SUITE_END() // Mapping
