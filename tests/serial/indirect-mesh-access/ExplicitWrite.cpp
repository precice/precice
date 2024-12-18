#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <numeric>
#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(IndirectMeshAccess)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

// Test case for a indirect mesh access on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. The test case here is the most basic variant in order
// use such a feature. SolverTwo defines the mesh whereas SolverOne writes
// indirectly to this mesh.
// nearest-neighbor-conservative-write
// Note that a data initialization is not possible with indirect and direct mesh
// access (see #1583)
BOOST_AUTO_TEST_CASE(ExplicitWrite)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 2;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0};
  std::vector<double>         writeData1({1, 2, -3, 4, 2});
  std::vector<double>         writeData2({-4, 12, 3, 5, 7});
  std::vector<double>         expectedData1({1, 4, -3, 4});
  std::vector<double>         expectedData2({-4, 19, 3, 5});

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

      // read data (not necessary here)
      // solve time step
      // write data:
      std::vector<double> tmpPositions = {0.0, -0.01, 0.01, 0.05, 0.1, 0.1, 0.1, 0.0, 0, 0.05};
      for (std::size_t i = 0; i < writeData1.size(); ++i) {
        std::vector<double> solverTwoCoord(dim);
        for (int d = 0; d < dim; ++d) {
          solverTwoCoord[d] = tmpPositions[i * dim + d];
        }
        if (time == 1) {
          couplingInterface.mapAndwriteData(otherMeshName, dataName, solverTwoCoord, {&writeData1[i], 1});
        } else if (time == 2) {
          couplingInterface.mapAndwriteData(otherMeshName, dataName, solverTwoCoord, {&writeData2[i], 1});
        } else {
          PRECICE_ASSERT(false);
        }
      }
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
    std::vector<double> readData(4, -1);

    double time = 0;
    // Initialize
    couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;
      couplingInterface.readData(meshName, dataName, ids, dt, readData);

      if (time == 1) {
        BOOST_TEST(std::accumulate(readData.begin(), readData.end(), 0.) == std::accumulate(writeData1.begin(), writeData1.end(), 0.));
        for (std::size_t i = 0; i < readData.size(); ++i) {
          BOOST_TEST(readData[i] == expectedData1[i]);
        }
      } else if (time == 2) {
        BOOST_TEST(std::accumulate(readData.begin(), readData.end(), 0.) == std::accumulate(writeData2.begin(), writeData2.end(), 0.));
        for (std::size_t i = 0; i < readData.size(); ++i) {
          BOOST_TEST(readData[i] == expectedData2[i]);
        }
      } else {
        PRECICE_ASSERT(false);
      }
      // solve time step
      // write data (not necessary here)
      couplingInterface.advance(dt);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
