#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test case for a direct mesh access by SolverTwo to a mesh defined
// by SolverOne. Both solvers read and write data to/from MeshOne.
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
BOOST_AUTO_TEST_CASE(DirectAccessReadWrite)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  if (context.isNamed("SolverOne")) {
    // Set up Participant
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim              = 2;
    const auto           providedMeshName = "MeshOne";
    const auto           readDataName     = "Forces";
    const auto           writeDataName    = "Velocities";
    BOOST_TEST(interface.getMeshDimensions(providedMeshName) == 2);

    std::vector<double> positions = std::vector<double>({0.5, 0.25});
    const int           meshSize  = positions.size() / dim;
    std::vector<int>    ids(meshSize, -1);
    interface.setMeshVertices(providedMeshName, positions, ids);

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Some dummy writeData
    std::vector<double> readData(ids.size(), -1);
    std::vector<double> writeData;
    for (int i = 0; i < meshSize; ++i)
      writeData.emplace_back(i + 5);

    int iterations = 0;
    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        // do nothing
      }

      interface.readData(providedMeshName, readDataName, ids, dt, readData);

      std::vector<double> expectedData = std::vector<double>({50});
      if (iterations == 0) {
        expectedData[0] = 0; // initial data
      }

      BOOST_TEST(precice::testing::equals(expectedData, readData));
      interface.writeData(providedMeshName, writeDataName, ids, writeData);
      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
      iterations++;
      if (interface.requiresReadingCheckpoint()) {
        // do nothing
      }
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim              = 2;
    const auto           receivedMeshName = "MeshOne";
    const auto           writeDataName    = "Forces";
    const auto           readDataName     = "Velocities";
    BOOST_TEST(interface.getMeshDimensions(receivedMeshName) == 2);

    std::array<double, dim * 2> boundingBox = std::array<double, dim * 2>{0.0, 1.0, 0.0, 1.0};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(receivedMeshName, boundingBox);

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int receivedMeshSize = interface.getMeshVertexSize(receivedMeshName);
    BOOST_TEST(receivedMeshSize == 1);

    // Allocate a vector containing the vertices
    std::vector<double> receivedMesh(receivedMeshSize * dim);
    std::vector<int>    receiveMeshIDs(receivedMeshSize, -1);
    interface.getMeshVertexIDsAndCoordinates(receivedMeshName, receiveMeshIDs, receivedMesh);

    // Allocate data to read and write
    std::vector<double> readData(receiveMeshIDs.size(), -1);
    std::vector<double> writeData;
    for (int i = 0; i < receivedMeshSize; ++i)
      writeData.emplace_back(i + 50);
    // Expected data = positions of the other participant's mesh
    const std::vector<double> expectedData = std::vector<double>({0.5, 0.25});
    BOOST_TEST(receivedMesh == expectedData);
    int iterations = 0;
    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        // do nothing
      }

      interface.readData(receivedMeshName, readDataName, receiveMeshIDs, dt, readData);

      std::vector<double> expectedData = std::vector<double>({5});
      if (iterations == 0) {
        expectedData[0] = 0; // initial data
      }

      BOOST_TEST(precice::testing::equals(expectedData, readData));
      interface.writeData(receivedMeshName, writeDataName, receiveMeshIDs, writeData);
      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
      iterations++;
      if (interface.requiresReadingCheckpoint()) {
        // do nothing
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // DirectMeshAccess

#endif // PRECICE_NO_MPI
