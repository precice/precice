#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
// Test case for a direct mesh access on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. As opposed to the 'boundingBoxExplicit' test case, this
// test case uses the same feature in an implicit setup.
BOOST_AUTO_TEST_CASE(Implicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  // Set up Solverinterface
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);
  constexpr int        dim = 2;

  if (context.isNamed("SolverOne")) {
    std::vector<double>         positions   = {0.1, 0.1, 0.2, 0.05, 0.1, 0.0, 0.3, 0.9};
    std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0};
    std::vector<int>            ownIDs(4, -1);

    auto ownMeshName   = "MeshOne";
    auto otherMeshName = "MeshTwo";
    auto ownDataName   = "Forces";
    auto otherDataName = "Velocities";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(ownMeshName) == 2);
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 2);

    // Define the own mesh
    couplingInterface.setMeshVertices(ownMeshName, positions, ownIDs);
    // TODO: Implement something in order to derive the bounding box from the mesh

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);

    couplingInterface.initialize();
    double dt = couplingInterface.getMaxTimeStepSize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int meshSize = couplingInterface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(meshSize == 3);

    // Allocate a vector containing the vertices
    std::vector<double> solverTwoMesh(meshSize * dim);
    std::vector<int>    otherIDs(meshSize);

    couplingInterface.getMeshVerticesAndIDs(otherMeshName, otherIDs, solverTwoMesh);
    // Some dummy writeData
    std::array<double, 3> writeData({1, 2, 3});

    // Expected data = positions of the other participant's mesh
    const std::vector<double> expectedData = {0.0, 0.0, 0.2, 0.3, 0.1, 0.1};
    BOOST_TEST(solverTwoMesh == expectedData);

    std::vector<double> readData(ownIDs.size(), -10);
    while (couplingInterface.isCouplingOngoing()) {
      if (couplingInterface.requiresWritingCheckpoint()) {
      }

      // Write data
      couplingInterface.writeData(otherMeshName, otherDataName, otherIDs, writeData);
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
      couplingInterface.readData(ownMeshName, ownDataName, ownIDs, dt, readData);
      if (couplingInterface.requiresReadingCheckpoint()) {
      }

      // Expected data according to the writeData
      std::vector<double> expectedData({10, 11, 12, 13});
      BOOST_TEST(precice::testing::equals(expectedData, readData));
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    std::vector<double>         positions = {0.0, 0.0, 0.2, 0.3, 0.1, 0.1};
    std::vector<int>            ownIDs(3, -1);
    std::array<double, dim * 2> boundingBox = {0.0, 2.0, 0.0, 2.0};

    // Query IDs
    auto ownMeshName   = "MeshTwo";
    auto otherMeshName = "MeshOne";
    auto ownDataName   = "Velocities";
    auto otherDataName = "Forces";

    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 2);
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(ownMeshName) == 2);

    // Define the mesh
    couplingInterface.setMeshVertices(ownMeshName, positions, ownIDs);
    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);
    // Initialize
    couplingInterface.initialize();
    double dt = couplingInterface.getMaxTimeStepSize();

    const int meshSize = couplingInterface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(meshSize == 4);

    // Allocate a vector containing the vertices
    std::vector<double> solverOneMesh(meshSize * dim);
    std::vector<int>    otherIDs(meshSize);

    couplingInterface.getMeshVerticesAndIDs(otherMeshName, otherIDs, solverOneMesh);
    // Some dummy writeData
    std::array<double, 4> writeData({10, 11, 12, 13});

    // Allocate data to read
    std::vector<double> readData(ownIDs.size(), -10);

    while (couplingInterface.isCouplingOngoing()) {
      if (couplingInterface.requiresWritingCheckpoint()) {
      }

      // Write data
      couplingInterface.writeData(otherMeshName, otherDataName, otherIDs, writeData);
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
      couplingInterface.readData(ownMeshName, ownDataName, ownIDs, dt, readData);
      if (couplingInterface.requiresReadingCheckpoint()) {
      }

      // Expected data according to the writeData
      std::vector<double> expectedData({1, 2, 3});
      BOOST_TEST(precice::testing::equals(expectedData, readData));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // AccessReceivedMesh

#endif // PRECICE_NO_MPI
