#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(IndirectMeshAccess)
// Test case for a indirect mesh access on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. The test case here is the most basic variant in order
// use such a feature. SolverTwo defines the mesh whereas SolverOne reads
// indirectly from this mesh.
BOOST_AUTO_TEST_CASE(ExplicitRead)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 2;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0};

  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    auto dataName      = "Velocities";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 2);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);

    couplingInterface.initialize();
    double dt = couplingInterface.getMaxTimeStepSize();

    // Allocate data to read
    while (couplingInterface.isCouplingOngoing()) {
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
      std::vector<double> expectedData({1, 2, 3, 4, 2, 4});
      std::vector<double> tmpPositions = {0.0, -0.01, 0.01, 0.05, 0.1, 0.1, 0.1, 0.0, 0, 0.05, 0.1, 0.0};

      for (std::size_t i = 0; i < expectedData.size(); ++i) {
        std::vector<double> solverTwoCoord(dim);
        double              value;
        for (int d = 0; d < dim; ++d) {
          solverTwoCoord[d] = tmpPositions[i * dim + d];
        }
        couplingInterface.readData(otherMeshName, dataName, solverTwoCoord, dt, {&value, 1});
        // Expected data according to the writeData
        BOOST_TEST(expectedData[i] == value);
      }
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName = "MeshTwo";
    auto dataName = "Velocities";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName));

    std::vector<double> positions = {0.0, 0.0, 0.0, 0.05, 0.1, 0.1, 0.1, 0.0};
    std::vector<int>    ids(4, -1);

    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // Some dummy readData
    std::array<double, 4> writeData({1, 2, 3, 4});

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

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
