#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(IndirectMeshAccess)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

// Test case for a indirect mesh access on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. The test case here is the most basic variant in order
// use such a feature. SolverTwo defines the mesh whereas SolverOne reads
// indirectly from this mesh.
// nearest-neighbor-consistent-read
BOOST_AUTO_TEST_CASE(ExplicitRead)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 2;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, -0.02, 1.0};

  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    auto dataName      = "Velocities";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 2);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);

    couplingInterface.initialize();

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // read data:
      std::vector<double> expectedData1({1, 2, 3, 4, 2, 4});
      std::vector<double> expectedData2({-10, -11, -12, -13, -11, -13});
      std::vector<double> tmpPositions = {0.0, -0.01, 0.01, 0.05, 0.1, 0.1, 0.1, 0.0, 0, 0.05, 0.1, 0.0};

      for (std::size_t i = 0; i < expectedData1.size(); ++i) {
        std::vector<double> solverTwoCoord(dim);
        double              value;
        for (int d = 0; d < dim; ++d) {
          solverTwoCoord[d] = tmpPositions[i * dim + d];
        }

        couplingInterface.mapAndReadData(otherMeshName, dataName, solverTwoCoord, dt, {&value, 1});
        // Expected data according to the writeData
        if (time == 1) {
          BOOST_TEST(expectedData1[i] == value);
        } else if (time == 2) {
          BOOST_TEST(expectedData2[i] == value);
        } else {
          PRECICE_ASSERT(false);
        }
      }

      // Check that we catch vertice not within the access region.
      {
        std::vector<double> solverTwoCoord(dim);
        double              value;
        for (int d = 0; d < dim; ++d) {
          solverTwoCoord[d] = 100 * tmpPositions[d];
        }
        // exceeding to the top
        BOOST_CHECK_THROW(couplingInterface.mapAndReadData(otherMeshName, dataName, solverTwoCoord, dt, {&value, 1}), ::precice::Error);

        // exceeding to the bottom
        for (int d = 0; d < dim; ++d) {
          solverTwoCoord[d] = -1000 * tmpPositions[d];
        }
        BOOST_CHECK_THROW(couplingInterface.mapAndReadData(otherMeshName, dataName, solverTwoCoord, dt, {&value, 1}), ::precice::Error);
      }

      // solve time step
      // write data (not necessary here)
      couplingInterface.advance(dt);
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
    std::array<double, 4> writeData1({1, 2, 3, 4});
    std::array<double, 4> writeData2({-10, -11, -12, -13});

    // Initialize
    couplingInterface.initialize();
    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;
      // read data (not necessary here)
      // solve time step
      // write data:
      if (time == 1) {
        couplingInterface.writeData(meshName, dataName, ids, writeData1);
      } else if (time == 2) {
        couplingInterface.writeData(meshName, dataName, ids, writeData2);
      } else {
        PRECICE_ASSERT(false);
      }
      couplingInterface.advance(dt);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
