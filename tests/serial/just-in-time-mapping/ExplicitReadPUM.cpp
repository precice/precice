#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <algorithm>
#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(JustInTimeMapping)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

// Test case for a just-in-time mapping on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. The test case here is the most basic variant using partition of unity to
// use such a feature. SolverTwo defines the mesh whereas SolverOne reads
// just-in-time from this mesh.
// pum-consistent-read with and without polynomial for scalar and vector data, reference is given by a conventional mapping
BOOST_AUTO_TEST_CASE(ExplicitReadPUM)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 2;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0};

  if (context.isNamed("SolverOne")) {
    auto meshName       = "MeshOne";
    auto otherMeshName  = "MeshTwo";
    auto meshBName      = "MeshOneB";
    auto otherMeshBName = "MeshTwoB";
    auto dataName       = "Temperature";
    auto vectorDataName = "Velocity";
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
    couplingInterface.setMeshVertices(meshBName, positions, ids);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);
    couplingInterface.setMeshAccessRegion(otherMeshBName, boundingBox);

    couplingInterface.initialize();

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // read data (for reference)
      std::vector<double> readData(meshSize);
      std::vector<double> justInTimeReadData(meshSize);
      std::vector<double> readVectorData(meshSize * dim);
      std::vector<double> justInTimeReadVectorData(meshSize * dim);

      // First, we check the separate polynomial PUM (scalar and vector)
      couplingInterface.readData(meshName, dataName, ids, dt, readData);
      couplingInterface.mapAndReadData(otherMeshName, dataName, positions, dt, justInTimeReadData);
      couplingInterface.readData(meshName, vectorDataName, ids, dt, readVectorData);
      couplingInterface.mapAndReadData(otherMeshName, vectorDataName, positions, dt, justInTimeReadVectorData);

      for (int r = 0; r < meshSize; ++r) {
        BOOST_TEST(readData[r] == justInTimeReadData[r]);
        BOOST_TEST(readVectorData[r * dim] == justInTimeReadVectorData[r * dim]);
        BOOST_TEST(readVectorData[r * dim + 1] == justInTimeReadVectorData[r * dim + 1]);
      }

      // Second, we check the no polynomial PUM (scalar and vector)
      couplingInterface.readData(meshBName, dataName, ids, dt, readData);
      couplingInterface.mapAndReadData(otherMeshBName, dataName, positions, dt, justInTimeReadData);
      couplingInterface.readData(meshBName, vectorDataName, ids, dt, readVectorData);
      couplingInterface.mapAndReadData(otherMeshBName, vectorDataName, positions, dt, justInTimeReadVectorData);

      for (int r = 0; r < meshSize; ++r) {
        BOOST_TEST(readData[r] == justInTimeReadData[r]);
        BOOST_TEST(readVectorData[r * dim] == justInTimeReadVectorData[r * dim]);
        BOOST_TEST(readVectorData[r * dim + 1] == justInTimeReadVectorData[r * dim + 1]);
      }

      // solve time step
      // write data (not necessary here)
      couplingInterface.advance(dt);
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName       = "MeshTwo";
    auto meshBName      = "MeshTwoB";
    auto dataName       = "Temperature";
    auto vectorDataName = "Velocity";
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
    couplingInterface.setMeshVertices(meshBName, positions, ids);
    // Some dummy readData
    std::vector<double> writeData1(100);
    std::vector<double> writeData2(100);
    std::vector<double> vectorData1(200);
    std::vector<double> vectorData2(200);

    // Linear filling for writeData1
    std::generate(writeData1.begin(), writeData1.end(), [n = 0]() mutable {
      return 0.5 * n++; // Simple increasing pattern, starting from 0
    });

    // Quadratic filling for writeData2
    std::generate(writeData2.begin(), writeData2.end(), [n = 0]() mutable {
      return -100.0 + 0.1 * (n * n++); // Quadratic pattern, n^2 scaled by 0.1
    });

    // vectorData1 with linearly increasing values
    for (int i = 0; i < 100; ++i) {
      vectorData1[2 * i]     = 0.5 * i; // x component
      vectorData1[2 * i + 1] = 0.2 * i; // y component
    }

    // vectorData2 with quadratic and linear patterns
    for (int i = 0; i < 100; ++i) {
      vectorData2[2 * i]     = 0.1 * i * i;   // x component quadratically increasing
      vectorData2[2 * i + 1] = 100 - 0.5 * i; // y component linearly decreasing
    }

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
        couplingInterface.writeData(meshName, vectorDataName, ids, vectorData1);
        couplingInterface.writeData(meshBName, dataName, ids, writeData1);
        couplingInterface.writeData(meshBName, vectorDataName, ids, vectorData1);
      } else if (time == 2) {
        couplingInterface.writeData(meshName, dataName, ids, writeData2);
        couplingInterface.writeData(meshName, vectorDataName, ids, vectorData2);
        couplingInterface.writeData(meshBName, dataName, ids, writeData2);
        couplingInterface.writeData(meshBName, vectorDataName, ids, vectorData2);
      } else {
        BOOST_TEST(false);
      }
      couplingInterface.advance(dt);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
