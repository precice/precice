#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

// Test case for a direct mesh access while also using waveform relaxation to test the integration of the two features.
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
BOOST_AUTO_TEST_CASE(DirectAccessWithWaveform)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  if (context.isNamed("SolverOne")) {
    // Set up Solverinterface
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    BOOST_TEST(interface.getDimensions() == 2);
    constexpr int dim         = 2;
    const int     ownMeshID   = interface.getMeshID("MeshOne");
    const int     otherMeshID = interface.getMeshID("MeshTwo");
    const int     readDataID  = interface.getDataID("Forces", ownMeshID);
    const int     writeDataID = interface.getDataID("Velocities", otherMeshID);

    std::vector<double> ownPositions = std::vector<double>({0.5, 0.25});
    std::vector<int>    ownIDs(ownPositions.size() / dim, -1);
    interface.setMeshVertices(ownMeshID, ownIDs.size(), ownPositions.data(), ownIDs.data());

    std::array<double, dim * 2> boundingBox = std::array<double, dim * 2>{0.0, 1.0, 0.0, 1.0};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(otherMeshID, boundingBox.data());

    double dt = interface.initialize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int otherMeshSize = interface.getMeshVertexSize(otherMeshID);
    BOOST_TEST(otherMeshSize == 1);

    std::vector<double> otherPositions(otherMeshSize * dim);
    std::vector<int>    otherIDs(otherMeshSize, -1);
    interface.getMeshVerticesAndIDs(otherMeshID, otherMeshSize, otherIDs.data(), otherPositions.data());

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

      interface.readBlockScalarData(readDataID, ownIDs.size(), ownIDs.data(), 0.5, readData.data());

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
      interface.writeBlockScalarData(writeDataID, otherIDs.size(), otherIDs.data(), writeData.data());
      dt = interface.advance(dt);
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
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    BOOST_TEST(interface.getDimensions() == 2);
    constexpr int dim         = 2;
    const int     meshID      = interface.getMeshID("MeshTwo");
    const int     writeDataID = interface.getDataID("Forces", meshID);
    const int     readDataID  = interface.getDataID("Velocities", meshID);

    std::vector<double> positions = std::vector<double>({0.5, 0.25});
    std::vector<int>    ids(positions.size() / dim, -1);
    interface.setMeshVertices(meshID, ids.size(), positions.data(), ids.data());

    double dt = interface.initialize();

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

      interface.readBlockScalarData(readDataID, ids.size(), ids.data(), 0.5, readData.data());

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
      interface.writeBlockScalarData(writeDataID, ids.size(), ids.data(), writeData.data());
      dt = interface.advance(dt);
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
