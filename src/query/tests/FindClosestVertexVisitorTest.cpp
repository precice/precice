#include <Eigen/Core>
#include <algorithm>
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "query/FindClosestVertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::query;

BOOST_AUTO_TEST_SUITE(QueryTests)

BOOST_AUTO_TEST_CASE(FindClosestVertexVisitor)
{
  PRECICE_TEST(1_rank);
  mesh::Mesh mesh("Mesh", 2, false, testing::nextMeshID());
  mesh.createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh.createVertex(Eigen::Vector2d(0.0, 5.0));
  FindClosestVertex find(Eigen::Vector2d(1, 0));
  bool              found = find(mesh);
  BOOST_TEST(found);
  mesh::Vertex &closestVertex = find.getClosestVertex();
  BOOST_TEST(testing::equals(closestVertex.getCoords(), Eigen::Vector2d(0, 0)));
  double distance = find.getEuclidianDistance();
  BOOST_TEST(distance == 1.0);
}

BOOST_AUTO_TEST_SUITE_END() // QueryTests
