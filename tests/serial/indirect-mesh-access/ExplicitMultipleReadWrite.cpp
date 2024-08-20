#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <numeric>
#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(IndirectMeshAccess)
// Test case for a indirect mesh access on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. SolverTwo defines the mesh whereas SolverOne reads
// indirectly from this mesh.
//
// This test case maps multiple data on the same mesh indirectly, ensuring that
// caching mechanisms for the mapping (which is shared across different data) works
// as expected.
//
// nearest-neighbor-consistent-read  (vector and scalar data)
// nearest-neighbor-conservative-write (vector and scalar data)
BOOST_AUTO_TEST_CASE(ExplicitMultipleReadWrite)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 3;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};

  // maps to 0, 1, 1 of the output mesh
  std::vector<double> tmpPositions1 = {0.0, 0.0, 0.01, 0.01, 0.01, 0.9, 0.0, 0.2, 1.0};
  std::vector<double> writeData1({27, 12, 8});
  std::vector<double> expectedData1({27, 20, 0, 0});

  // std::vector<double> tmpPositions2 = {0.0, -0.01, 0.01, 0.05, 0.1, 0.1, 0.1, 0.0, 0, 0.05};
  std::vector<double> writeData2({27, 12, 8});
  std::vector<double> expectedData2({27, 20, 0, 0});

  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    auto dataName      = "Velocities";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 3);

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
      for (std::size_t i = 0; i < writeData1.size(); ++i) {
        std::vector<double> solverTwoCoord(dim);
        if (time == 1) {
          for (int d = 0; d < dim; ++d) {
            solverTwoCoord[d] = tmpPositions1[i * dim + d];
          }
          couplingInterface.mapAndwriteData(otherMeshName, dataName, solverTwoCoord, {&writeData1[i], 1});
        } else if (time == 2) {
          for (int d = 0; d < dim; ++d) {
            solverTwoCoord[d] = tmpPositions1[i * dim + d];
          }
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
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName) == 3);

    std::vector<double> positions = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};
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
