#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/impl/CreateClustering.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;
using namespace precice::mapping;
using namespace precice::testing;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(PartitionOfUnityClustering)

BOOST_AUTO_TEST_SUITE(Serial)

BOOST_AUTO_TEST_CASE(createClustering2D)
{
  PRECICE_TEST(1_rank);

  int meshDimension = 2;
  // Generate the meshes
  mesh::PtrMesh inMesh  = std::make_shared<Mesh>("inMesh", meshDimension, testing::nextMeshID());
  mesh::PtrMesh outMesh = std::make_shared<Mesh>("outMesh", meshDimension, testing::nextMeshID());

  // Create triangular shaped matching meshes
  for (unsigned int i = 0; i < 10; ++i) {
    for (unsigned int j = i; j < 10; ++j) {
      inMesh->createVertex(Eigen::Vector2d(static_cast<double>(i), static_cast<double>(j)));
      outMesh->createVertex(Eigen::Vector2d(static_cast<double>(i), static_cast<double>(j)));
    }
  }
  double       relativeOverlap      = 0.3;
  unsigned int verticesPerPartition = 10;
  bool         projectToInput       = false;
  {
    auto [averagePartitionRadius, centerCandidates] = impl::createClustering(inMesh, outMesh, relativeOverlap, verticesPerPartition, projectToInput);
    BOOST_TEST(averagePartitionRadius == 2.2360679774997898);
    BOOST_TEST(centerCandidates.size() == 19);
  }

  {
    projectToInput                                  = true;
    auto [averagePartitionRadius, centerCandidates] = impl::createClustering(inMesh, outMesh, relativeOverlap, verticesPerPartition, projectToInput);
    BOOST_TEST(averagePartitionRadius == 2.2360679774997898);
    BOOST_TEST(centerCandidates.size() == 25);
  }
}

BOOST_AUTO_TEST_CASE(createClustering3D)
{
  PRECICE_TEST(1_rank);

  int meshDimension = 3;
  // Generate the meshes
  mesh::PtrMesh inMesh  = std::make_shared<Mesh>("inMesh", meshDimension, testing::nextMeshID());
  mesh::PtrMesh outMesh = std::make_shared<Mesh>("outMesh", meshDimension, testing::nextMeshID());

  // Create triangular shaped matching meshes
  for (unsigned int i = 0; i < 10; ++i) {
    for (unsigned int j = i; j < 10; ++j) {
      for (unsigned int k = i; k < 10; ++k) {
        inMesh->createVertex(Eigen::Vector3d(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)));
        outMesh->createVertex(Eigen::Vector3d(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)));
      }
    }
  }
  double       relativeOverlap      = 0.3;
  unsigned int verticesPerPartition = 10;
  bool         projectToInput       = false;
  {
    auto [averagePartitionRadius, centerCandidates] = impl::createClustering(inMesh, outMesh, relativeOverlap, verticesPerPartition, projectToInput);
    BOOST_TEST(averagePartitionRadius == 1.4142135623730951);
    BOOST_TEST(centerCandidates.size() == 188);
  }

  {
    projectToInput                                  = true;
    auto [averagePartitionRadius, centerCandidates] = impl::createClustering(inMesh, outMesh, relativeOverlap, verticesPerPartition, projectToInput);
    BOOST_TEST(averagePartitionRadius == 1.4142135623730951);
    BOOST_TEST(centerCandidates.size() == 222);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
