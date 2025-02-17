#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(JustInTimeMapping)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

// Test case for a just-in-time mapping on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. The read direction is tested by comparing the just-in-time mapping against
// the conventional preCICE functions
// implicit coupling, nearest-neighbor, just-in-time mapping (write-conservative, read-consistent)
BOOST_AUTO_TEST_CASE(Implicit)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);
  constexpr int        dim = 2;

  // Allocate read and write data
  std::vector<double> writeDataForce({2.2, 3.3, 4.1});

  if (context.isNamed("SolverOne")) {
    std::array<double, dim * 2> boundingBox        = {0.0, 1.0, 0.0, 1.0};
    std::vector<double>         tmpMeshCoordsWrite = {0.1, 0.1, 0.2, 0.05};
    std::vector<double>         tmpMeshCoordsRead  = {0.2, 0.3};
    std::vector<int>            testIDs            = {-1};
    std::vector<int>            test2IDs           = {-1, -1};
    std::vector<double>         testReadData({-1.});

    auto otherMeshName = "MeshTwo";
    auto testMeshName  = "TestMesh";
    auto testMesh2Name = "TestMesh2";
    auto ownDataName   = "Velocities";
    auto otherDataName = "Forces";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == dim);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);
    couplingInterface.setMeshVertices(testMeshName, tmpMeshCoordsRead, testIDs);
    couplingInterface.setMeshVertices(testMesh2Name, tmpMeshCoordsWrite, test2IDs);

    couplingInterface.initialize();

    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    // more relevant for direct access, but should work
    const int meshSize = couplingInterface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(meshSize == 3);

    // Some dummy writeData
    std::vector<double> writeData({4, 2});
    std::vector<double> readData(tmpMeshCoordsRead.size() / dim, -10);

    int    iteration   = 0;
    double timestep    = 0;
    double oldTimestep = 0;
    while (couplingInterface.isCouplingOngoing()) {
      if (couplingInterface.requiresWritingCheckpoint()) {
        oldTimestep = timestep;
      }
      double dt = couplingInterface.getMaxTimeStepSize();
      ++iteration;
      ++timestep;
      // read data
      couplingInterface.mapAndReadData(otherMeshName, otherDataName, tmpMeshCoordsRead, dt, readData);
      // TODO: prevent ID access
      couplingInterface.readData(testMeshName, otherDataName, testIDs, dt, testReadData);
      BOOST_TEST(precice::testing::equals(testReadData, readData));

      // solve timestep
      // write data
      std::transform(writeData.begin(), writeData.end(), writeData.begin(), [&](auto &w) { return timestep * 500 + 50 * std::pow(10, -iteration); });
      couplingInterface.writeAndMapData(otherMeshName, ownDataName, tmpMeshCoordsWrite, writeData);
      couplingInterface.writeData(testMesh2Name, ownDataName, test2IDs, writeData);
      couplingInterface.advance(dt);

      if (couplingInterface.requiresReadingCheckpoint()) {
        timestep = oldTimestep;
      } else {
        iteration = 0;
      }
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    std::vector<double> positions = {0.0, 0.0, 0.2, 0.3, 0.1, 0.1};
    std::vector<int>    ownIDs(3, -1);
    std::vector<int>    testIDs(3, -1);

    // Query IDs
    auto ownMeshName   = "MeshTwo";
    auto testMeshName  = "MeshTwoDuplicate";
    auto ownDataName   = "Forces";
    auto otherDataName = "Velocities";

    BOOST_REQUIRE(couplingInterface.getMeshDimensions(ownMeshName) == dim);

    // Define the mesh
    couplingInterface.setMeshVertices(ownMeshName, positions, ownIDs);
    couplingInterface.setMeshVertices(testMeshName, positions, testIDs);

    // Initialize
    couplingInterface.initialize();

    std::vector<double> readData(ownIDs.size(), -10);
    std::vector<double> testReadData(testIDs.size(), -10);

    int    iteration   = 0;
    double timestep    = 0;
    double oldTimestep = 0;
    while (couplingInterface.isCouplingOngoing()) {
      if (couplingInterface.requiresWritingCheckpoint()) {
        oldTimestep = timestep;
      }
      double dt = couplingInterface.getMaxTimeStepSize();
      ++iteration;
      ++timestep;

      // read data
      couplingInterface.readData(ownMeshName, otherDataName, ownIDs, dt, readData);
      couplingInterface.readData(testMeshName, otherDataName, testIDs, dt, testReadData);
      BOOST_TEST(precice::testing::equals(testReadData, readData));
      // solve timestep
      // Write data
      std::transform(writeDataForce.begin(), writeDataForce.end(), writeDataForce.begin(), [&](auto &w) { return timestep * 5 + std::pow(10, -iteration); });
      couplingInterface.writeData(ownMeshName, ownDataName, ownIDs, writeDataForce);
      couplingInterface.advance(dt);

      if (couplingInterface.requiresReadingCheckpoint()) {
        timestep = oldTimestep;
      } else {
        iteration = 0;
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
