#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
// Test case for a direct mesh access on one participant to a mesh defined
// by another participant. Here, we check that preCICE correctly throws an
// error in case API functions are used, but the flag was not set in the config
// We also check that getMeshVertexSize can still be called and returns the
// amount of vertices in the remote mesh
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ExplicitReadNoFlag)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  std::vector<double> positions = {0.0, 0.0, 0.0, 0.05, 0.1, 0.1, 0.1, 0.0, 5.0, 5.0};
  std::vector<int>    ids(5, -1);

  constexpr int               dim         = 2;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0};

  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 2);

    // Check for thwrowing when trying to set the region
    BOOST_CHECK_THROW(couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox), ::precice::Error);

    couplingInterface.initialize();
    double dt = couplingInterface.getMaxTimeStepSize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int meshSize = couplingInterface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(meshSize == (ids.size()));

    // Allocate a vector containing the vertices
    std::vector<double> solverTwoMesh(meshSize * dim);
    BOOST_CHECK_THROW(couplingInterface.getMeshVertexIDsAndCoordinates(otherMeshName, ids, solverTwoMesh), ::precice::Error);
    while (couplingInterface.isCouplingOngoing()) {
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName = "MeshTwo";
    auto dataName = "Velocities";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName));

    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // Some dummy readData
    std::array<double, 5> writeData({1, 2, 3, 4, 5});

    // Initialize
    couplingInterface.initialize();
    double dt = couplingInterface.getMaxTimeStepSize();
    while (couplingInterface.isCouplingOngoing()) {

      couplingInterface.writeData(meshName, dataName, ids, writeData);
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // AccessReceivedMesh

#endif // PRECICE_NO_MPI
