#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test case for a direct mesh access while also using waveform relaxation to test the integration of the two features.
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
BOOST_AUTO_TEST_CASE(DirectAccessWithWaveform)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  if (context.isNamed("SolverOne")) {
    // Set up Participant
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim           = 2;
    const auto           ownMeshName   = "MeshOne";
    const auto           otherMeshName = "MeshTwo";
    const auto           readDataName  = "Forces";
    const auto           writeDataName = "Velocities";
    BOOST_REQUIRE(interface.getMeshDimensions(ownMeshName) == 2);
    BOOST_REQUIRE(interface.getMeshDimensions(otherMeshName) == 2);

    std::vector<double> ownPositions = std::vector<double>({0.5, 0.25});
    std::vector<int>    ownIDs(ownPositions.size() / dim, -1);
    interface.setMeshVertices(ownMeshName, ownPositions, ownIDs);

    std::array<double, dim * 2> boundingBox = std::array<double, dim * 2>{0.0, 1.0, 0.0, 1.0};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(otherMeshName, boundingBox);

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int otherMeshSize = interface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(otherMeshSize == 1);

    std::vector<double> otherPositions(otherMeshSize * dim);
    std::vector<int>    otherIDs(otherMeshSize, -1);
    interface.getMeshVerticesAndIDs(otherMeshName, otherIDs, otherPositions);

    // Some dummy writeData
    std::vector<double> readData(ownIDs.size(), -1);
    std::vector<double> writeData;
    for (int i = 0; i < otherMeshSize; ++i) {
      writeData.emplace_back(5);
    }

    int iterations = 0;
    int timeWindow = 0;

    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        // do nothing
      }

      interface.readData(ownMeshName, readDataName, ownIDs, 0.5, readData);

      std::vector<double> expectedData = std::vector<double>({-1});

      if (timeWindow == 0) {
        if (iterations == 0) {
          expectedData[0] = 0; // initial data
        } else {
          expectedData[0] = 25;
        }
      } else if (timeWindow == 1) {
        if (iterations == 0) {
          expectedData[0] = 50; // constant initial guess
        } else {
          expectedData[0] = 55;
        }
      }

      BOOST_TEST(precice::testing::equals(expectedData, readData));
      interface.writeData(otherMeshName, writeDataName, otherIDs, writeData);
      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
      iterations++;
      if (interface.requiresReadingCheckpoint()) {
        // do nothing
      }
      if (interface.isTimeWindowComplete()) {
        timeWindow++;
        iterations = 0;
        for (int i = 0; i < otherMeshSize; ++i) {
          writeData[i] = writeData[i] + 2;
        }
      }
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim         = 2;
    const auto           meshID      = "MeshTwo";
    const auto           writeDataID = "Forces";
    const auto           readDataID  = "Velocities";
    BOOST_REQUIRE(interface.getMeshDimensions(meshID) == 2);

    std::vector<double> positions = std::vector<double>({0.5, 0.25});
    std::vector<int>    ids(positions.size() / dim, -1);
    interface.setMeshVertices(meshID, positions, ids);

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Allocate data to read and write
    std::vector<double> readData(ids.size(), -1);
    std::vector<double> writeData;
    for (unsigned int i = 0; i < ids.size(); ++i) {
      writeData.emplace_back(50);
    }

    int iterations = 0;
    int timeWindow = 0;

    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        // do nothing
      }

      interface.readData(meshID, readDataID, ids, 0.5, readData);

      std::vector<double> expectedData = std::vector<double>({-1});

      if (timeWindow == 0) {
        if (iterations == 0) {
          expectedData[0] = 0; // initial data
        } else {
          expectedData[0] = 2.5;
        }
      } else if (timeWindow == 1) {
        if (iterations == 0) {
          expectedData[0] = 5; // constant initial guess
        } else {
          expectedData[0] = 6;
        }
      }

      BOOST_TEST(precice::testing::equals(expectedData, readData));
      interface.writeData(meshID, writeDataID, ids, writeData);
      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
      iterations++;
      if (interface.requiresReadingCheckpoint()) {
        // do nothing
      }
      if (interface.isTimeWindowComplete()) {
        timeWindow++;
        iterations = 0;
        for (std::size_t i = 0; i < ids.size(); ++i) {
          writeData[i] = writeData[i] + 10;
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // DirectMeshAccess

#endif // PRECICE_NO_MPI
