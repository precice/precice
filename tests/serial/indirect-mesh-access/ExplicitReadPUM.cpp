#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <algorithm>
#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(IndirectMeshAccess)
// Test case for a indirect mesh access on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. The test case here is the most basic variant using partition of unity to
// use such a feature. SolverTwo defines the mesh whereas SolverOne reads
// indirectly from this mesh.
// pum-consistent-read
BOOST_AUTO_TEST_CASE(ExplicitReadPUM)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 2;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0};

  if (context.isNamed("SolverOne")) {
    auto meshName      = "MeshOne";
    auto otherMeshName = "MeshTwo";
    auto dataName      = "Velocities";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 2);

    int                 meshSize = 10;
    std::vector<int>    ids(meshSize);
    std::vector<double> positions;
    positions.reserve(dim * meshSize);

    // Adding 10 test positions
    double offset = 0.05;         // small offset to ensure no overlap with grid points
    for (int k = 0; k < 5; ++k) { // Two loops, but only 5 iterations each to add 10 new positions
      for (int l = 0; l < 2; ++l) {
        positions.emplace_back(k * 0.2 + offset); // Increment x by 0.2 each time, starting from 0.1
        positions.emplace_back(l * 0.8 + offset); // Two y positions: near 0.1 and near 0.9
      }
    }
    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshVertices(meshName, positions, ids);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);

    couplingInterface.initialize();

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // read data (for reference)
      std::vector<double> readData(meshSize);
      std::vector<double> indirectReadData(meshSize);
      couplingInterface.readData(meshName, dataName, ids, dt, readData);
      couplingInterface.mapAndreadData(otherMeshName, dataName, positions, dt, indirectReadData);

      for(int r = 0; r<meshSize;++r)
      BOOST_TEST(readData[r] == indirectReadData[r]);
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

    std::vector<double> positions;
    positions.reserve(200);
    for (int i = 0; i < 10; ++i) {
      for (int j = 0; j < 10; ++j) {
        positions.emplace_back(i * 0.1);
        positions.emplace_back(j * 0.1);
      }
    }
    std::vector<int> ids(100, -1);
    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // Some dummy readData
    std::vector<double> writeData1(100);
    std::vector<double> writeData2(100);

    // Linear filling for writeData1
    std::generate(writeData1.begin(), writeData1.end(), [n = 0]() mutable {
      return 0.5 * n++; // Simple increasing pattern, starting from 0
    });

    // Quadratic filling for writeData2
    std::generate(writeData2.begin(), writeData2.end(), [n = 0]() mutable {
      return -100.0 + 0.1 * (n * n++); // Quadratic pattern, n^2 scaled by 0.1
    });

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
