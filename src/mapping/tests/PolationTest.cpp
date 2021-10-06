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
  mesh::Vertex vertex(Eigen::Vector3d(1.0, 2.0, 0.0), 0);

  Polation polation(vertex);

  std::vector<int>    expectedIndices = {0};
  std::vector<double> expectedWeights = {1.0};

  BOOST_TEST(polation.getWeightedElements().size() == 1);
  BOOST_TEST(polation.isInterpolation());

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
  mesh::Edge   edge(v1, v2, 0);

  Eigen::Vector3d location(0.0, 0.4, 0.0);

  Polation polation(location, edge);

  std::vector<int>    expectedIndices = {0, 1};
  std::vector<double> expectedWeights = {0.8, 0.2};

  BOOST_TEST(polation.getWeightedElements().size() == 2);
  BOOST_TEST(polation.isInterpolation());

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
  mesh::Edge     e1(v1, v2, 0);
  mesh::Edge     e2(v2, v3, 1);
  mesh::Edge     e3(v1, v3, 2);
  mesh::Triangle triangle(e1, e2, e3, 0);
  triangle.computeNormal();

  Eigen::Vector3d location(1.0, 0.6, 0.0);

  Polation polation(location, triangle);

  std::vector<int>    expectedIndices = {0, 1, 2};
  std::vector<double> expectedWeights = {0.35, 0.35, 0.3};

  BOOST_TEST(polation.getWeightedElements().size() == 3);
  BOOST_TEST(polation.isInterpolation());

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
  mesh::Edge   edge(v1, v2, 0);

  Eigen::Vector3d location(0.0, 3.0, 0.0);

  Polation polation(location, edge);

  std::vector<int>    expectedIndices = {0, 1};
  std::vector<double> expectedWeights = {-0.5, 1.5};

  BOOST_TEST(polation.getWeightedElements().size() == 2);
  BOOST_TEST(not polation.isInterpolation());

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
  mesh::Edge     e1(v1, v2, 0);
  mesh::Edge     e2(v2, v3, 1);
  mesh::Edge     e3(v1, v3, 2);
  mesh::Triangle triangle(e1, e2, e3, 0);
  triangle.computeNormal();

  Eigen::Vector3d location(4.0, 0.6, 0.0);

  Polation polation(location, triangle);

  std::vector<int>    expectedIndices = {0, 1, 2};
  std::vector<double> expectedWeights = {-1.15, 1.85, 0.3};

  BOOST_TEST(polation.getWeightedElements().size() == 3);
  BOOST_TEST(not polation.isInterpolation());

  for (size_t i = 0; i < polation.getWeightedElements().size(); ++i) {

    BOOST_TEST(polation.getWeightedElements().at(i).weight == expectedWeights.at(i));
    BOOST_TEST(polation.getWeightedElements().at(i).vertexID == expectedIndices.at(i));
  }
}

BOOST_AUTO_TEST_SUITE_END() // Interpolation
BOOST_AUTO_TEST_SUITE_END() // Mapping