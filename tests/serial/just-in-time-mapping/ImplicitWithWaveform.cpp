#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test case for a just-in-time mapping while also using waveform relaxation to test the integration of the two features.
// Here, we test the substeps = true
// The test case uses the "conventional" preCICE API, i.e., usual mappings as a reference solution
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(JustInTimeMapping)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ImplicitWithWaveform)
{
  PRECICE_TEST();
  constexpr int dim = 3;

  if (context.isNamed("SolverOne")) {
    // Set up Participant
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    const auto           ownMeshName    = "MeshOne";
    const auto           otherMeshName  = "MeshTwo";
    const auto           readDataName   = "Forces";
    const auto           writeDataName  = "Velocities";
    const auto           readDataName1  = "O2";
    const auto           writeDataName1 = "H2";
    BOOST_REQUIRE(interface.getMeshDimensions(ownMeshName) == 3);
    BOOST_REQUIRE(interface.getMeshDimensions(otherMeshName) == 3);

    std::vector<double> positions;
    double              offset = 0.03; // small offset to ensure no overlap with grid points
    for (int k = 0; k < 5; ++k) {      // Three loops, but only 5 iterations each to add 10 new positions
      for (int l = 0; l < 2; ++l) {
        for (int m = 0; m < 3; ++m) {
          positions.emplace_back(k * 0.2 + offset);  // Increment x by 0.2 each time, starting from 0.1
          positions.emplace_back(l * 0.8 + offset);  // Two y positions: near 0.1 and near 0.9
          positions.emplace_back(m * 0.15 + offset); // Three z positions
        }
      }
    }
    std::size_t      size = positions.size() / dim;
    std::vector<int> ownIDs(size, -1);
    interface.setMeshVertices(ownMeshName, positions, ownIDs);

    std::array<double, dim * 2> boundingBox = std::array<double, dim * 2>{0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(otherMeshName, boundingBox);

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Some dummy writeData
    std::vector<double> readData(size, -1);
    std::vector<double> referenceReadData(size, -1);
    std::vector<double> readDataTime(size, -1);
    std::vector<double> referenceReadDataTime(size, -1);
    std::vector<double> writeData(size, -1);

    std::vector<double> readData1(size * dim, -1);
    std::vector<double> referenceReadData1(size * dim, -1);
    std::vector<double> writeData1(size * dim, -1);

    int iterations = 0;
    int timeWindow = 0;
    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        // do nothing
      }

      // We read here alternating in time to check the cache
      interface.mapAndReadData(otherMeshName, readDataName, {positions.data(), dim * (size / 2)}, 0.1 * dt, {readData.data(), size / 2});
      interface.readData(ownMeshName, readDataName, ownIDs, 0.1 * dt, referenceReadData);
      interface.mapAndReadData(otherMeshName, readDataName, positions, 0.15 * dt, readDataTime);
      interface.readData(ownMeshName, readDataName, ownIDs, 0.15 * dt, referenceReadDataTime);
      interface.mapAndReadData(otherMeshName, readDataName, {positions.data() + dim * (size / 2), dim * (size / 2)}, 0.1 * dt, {readData.data() + size / 2, size / 2});

      BOOST_TEST(referenceReadData == readData, boost::test_tools::per_element());
      BOOST_TEST(referenceReadDataTime == readDataTime, boost::test_tools::per_element());

      interface.mapAndReadData(otherMeshName, readDataName1, positions, 0.1 * dt, readData1);
      interface.readData(ownMeshName, readDataName1, ownIDs, 0.1 * dt, referenceReadData1);

      BOOST_TEST(referenceReadData1 == readData1, boost::test_tools::per_element());
      // Some quadratic filling including contributions from the iteration and timeWindow
      std::generate(writeData.begin(), writeData.end(), [n = 0, timeWindow, iterations]() mutable {
        return 28 * n * n++ + timeWindow - iterations * n;
      });

      std::generate(writeData1.begin(), writeData1.end(), [n = 0, timeWindow, iterations]() mutable {
        return 4 * n++ + timeWindow * n - iterations;
      });
      // Just in time variant
      interface.writeAndMapData(otherMeshName, writeDataName, positions, writeData);
      interface.writeData(ownMeshName, writeDataName, ownIDs, writeData);
      interface.writeAndMapData(otherMeshName, writeDataName1, positions, writeData1);
      interface.writeData(ownMeshName, writeDataName1, ownIDs, writeData1);
      interface.advance(0.2 * dt);
      iterations++;
      if (interface.requiresReadingCheckpoint()) {
        // do nothing
      }
      if (interface.isTimeWindowComplete()) {
        timeWindow++;
        iterations = 0;
        dt         = interface.getMaxTimeStepSize();
      }
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    const auto           meshID       = "MeshTwo";
    const auto           refMesh      = "MeshThree";
    const auto           writeDataID  = "Forces";
    const auto           readDataID   = "Velocities";
    const auto           writeDataID1 = "O2";
    const auto           readDataID1  = "H2";

    std::vector<double> positions;
    int                 size = 125;
    positions.reserve(size * dim);
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
        for (int k = 0; k < 5; ++k) {
          positions.emplace_back(i * 0.2);
          positions.emplace_back(j * 0.2);
          positions.emplace_back(k * 0.2);
        }
      }
    }
    std::vector<int> ids(positions.size() / dim, -1);
    std::vector<int> refIDs(positions.size() / dim, -1);
    interface.setMeshVertices(meshID, positions, ids);
    interface.setMeshVertices(refMesh, positions, refIDs);

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Allocate data to read and write
    std::vector<double> readData(ids.size(), -1);
    std::vector<double> refReadData(ids.size(), -1);
    std::vector<double> writeData(ids.size(), -1);
    std::vector<double> readData1(ids.size() * dim, -1);
    std::vector<double> refReadData1(ids.size() * dim, -1);
    std::vector<double> writeData1(ids.size() * dim, -1);

    int iterations = 0;
    int timeWindow = 0;

    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        // do nothing
      }
      interface.readData(meshID, readDataID, ids, 0.25 * dt, readData);
      interface.readData(refMesh, readDataID, refIDs, 0.25 * dt, refReadData);
      interface.readData(meshID, readDataID1, ids, 0.25 * dt, readData1);
      interface.readData(refMesh, readDataID1, refIDs, 0.25 * dt, refReadData1);
      BOOST_TEST(refReadData == readData, boost::test_tools::per_element());
      BOOST_TEST(refReadData1 == readData1, boost::test_tools::per_element());

      // Some quadratic filling including contributions from the iteration and timeWindow
      std::generate(writeData.begin(), writeData.end(), [n = 0, timeWindow, iterations]() mutable {
        return 0.8 * n * n++ + timeWindow - iterations * n;
      });
      std::generate(writeData1.begin(), writeData1.end(), [n = 0, timeWindow, iterations]() mutable {
        return 1e-3 * n * n * n++ + timeWindow - iterations * n;
      });

      interface.writeData(meshID, writeDataID, ids, writeData);
      interface.writeData(meshID, writeDataID1, ids, writeData1);
      interface.advance(0.5 * dt);
      iterations++;
      if (interface.requiresReadingCheckpoint()) {
        // do nothing
      }
      if (interface.isTimeWindowComplete()) {
        timeWindow++;
        iterations = 0;
        dt         = interface.getMaxTimeStepSize();
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Just-in-time mapping

#endif // PRECICE_NO_MPI
