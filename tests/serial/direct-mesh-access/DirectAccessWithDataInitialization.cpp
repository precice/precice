#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test case for a direct mesh access while also using data initialization to test the integration of the two features.
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
BOOST_AUTO_TEST_CASE(DirectAccessWithDataInitialization)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  if (context.isNamed("SolverOne")) {
    // Set up Participant
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim         = 2;
    const auto           ownMeshID   = "MeshOne";
    const auto           otherMeshID = "MeshTwo";
    const auto           readDataID  = "Forces";
    const auto           writeDataID = "Velocities";
    BOOST_REQUIRE(interface.getMeshDimensions(ownMeshID) == 2);
    BOOST_REQUIRE(interface.getMeshDimensions(otherMeshID) == 2);

    std::vector<double> ownPositions = std::vector<double>({0.5, 0.25});
    std::vector<int>    ownIDs(ownPositions.size() / dim, -1);
    interface.setMeshVertices(ownMeshID, ownPositions, ownIDs);

    std::array<double, dim * 2> boundingBox = std::array<double, dim * 2>{0.0, 1.0, 0.0, 1.0};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(otherMeshID, boundingBox);

    std::vector<double> readData(ownIDs.size(), -1);
    int                 otherMeshSize = 1; // @todo hard-coded, because we cannot read this from preCICE before interface.initialize(). See https://github.com/precice/precice/issues/1583.
    std::vector<double> writeData(otherMeshSize, -1);

    BOOST_TEST(!interface.requiresInitialData());

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    BOOST_TEST(otherMeshSize == interface.getMeshVertexSize(otherMeshID)); // @todo would need to know this already earlier (see above).
    BOOST_TEST(otherMeshSize == 1);

    std::vector<double> otherPositions(otherMeshSize * dim);
    std::vector<int>    otherIDs(otherMeshSize, -1);
    interface.getMeshVerticesAndIDs(otherMeshID, otherIDs, otherPositions);

    // writeData for first window
    for (int i = 0; i < otherMeshSize; ++i) {
      writeData[i] = 5;
    }

    int iterations = 0;
    int timeWindow = 0;

    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        // do nothing
      }

      interface.readData(ownMeshID, readDataID, ownIDs, dt, readData);

      std::vector<double> expectedData = std::vector<double>({-1});

      if (timeWindow == 0) {
        if (iterations == 0) {
          expectedData[0] = 20; // initial data
        } else {
          expectedData[0] = 50;
        }
      } else if (timeWindow == 1) {
        if (iterations == 0) {
          expectedData[0] = 50; // constant initial guess
        } else {
          expectedData[0] = 60;
        }
      }

      BOOST_TEST(precice::testing::equals(expectedData, readData));
      interface.writeData(otherMeshID, writeDataID, otherIDs, writeData);
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
    constexpr int        dim           = 2;
    const auto           meshName      = "MeshTwo";
    const auto           writeDataName = "Forces";
    const auto           readDataName  = "Velocities";
    BOOST_REQUIRE(interface.getMeshDimensions(meshName) == 2);

    std::vector<double> positions = std::vector<double>({0.5, 0.25});
    std::vector<int>    ids(positions.size() / dim, -1);
    interface.setMeshVertices(meshName, positions, ids);

    std::vector<double> readData(ids.size(), -1);
    std::vector<double> writeData;

    // writeData for initialization
    for (std::size_t i = 0; i < ids.size(); ++i) {
      writeData.emplace_back(20);
    }

    if (interface.requiresInitialData()) {
      interface.writeData(meshName, writeDataName, ids, writeData);
    }

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // writeData for first window
    for (unsigned int i = 0; i < ids.size(); ++i) {
      writeData[i] = 50;
    }

    int iterations = 0;
    int timeWindow = 0;

    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        // do nothing
      }

      interface.readData(meshName, readDataName, ids, dt, readData);

      std::vector<double> expectedData = std::vector<double>({-1});

      if (timeWindow == 0) {
        if (iterations == 0) {
          // expectedData[0] = 2; // initial data, @todo Currently not possible to provide initial data != 0 via direct access. See https://github.com/precice/precice/issues/1583
          expectedData[0] = 0; // initial data
        } else {
          expectedData[0] = 5;
        }
      } else if (timeWindow == 1) {
        if (iterations == 0) {
          expectedData[0] = 5; // initial guess
        } else {
          expectedData[0] = 7;
        }
      }

      BOOST_TEST(precice::testing::equals(expectedData, readData));
      interface.writeData(meshName, writeDataName, ids, writeData);
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
