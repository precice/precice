#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <numeric>
#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(JustInTimeMapping)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

// Test case for a just-in-time mapping on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. The test case here is the most basic variant in order
// use such a feature. SolverTwo defines the mesh whereas SolverOne writes
// just-in-time to this mesh.
// pum-rbf-conservative-write
// Note that a data initialization is not possible with just-in-time mapping and direct mesh
// access (see #1583)
BOOST_AUTO_TEST_CASE(ExplicitWritePUM)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 3;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};

  if (context.isNamed("SolverOne")) {
    auto meshName       = "MeshOne";
    auto otherMeshName  = "MeshTwo";
    auto otherMeshBName = "MeshTwoB";
    auto dataName       = "Temperature";
    auto vectorDataName = "Velocity";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 3);

    int                 meshSize = 30;
    std::vector<int>    ids(meshSize);
    std::vector<double> positions;
    positions.reserve(dim * meshSize);

    // Adding 30 test positions
    double offset = 0.05;         // small offset to ensure no overlap with grid points
    for (int k = 0; k < 5; ++k) { // Three loops, but only 5 iterations each to add 10 new positions
      for (int l = 0; l < 2; ++l) {
        for (int m = 0; m < 3; ++m) {
          positions.emplace_back(k * 0.2 + offset);  // Increment x by 0.2 each time, starting from 0.1
          positions.emplace_back(l * 0.8 + offset);  // Two y positions: near 0.1 and near 0.9
          positions.emplace_back(m * 0.15 + offset); // Three z positions
        }
      }
    }
    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // couplingInterface.setMeshVertices(meshBName, positions, ids);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);
    couplingInterface.setMeshAccessRegion(otherMeshBName, boundingBox);

    couplingInterface.initialize();

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // write data (for reference)
      std::vector<double> writeData(meshSize);
      std::vector<double> justInTimeWriteData(meshSize);
      std::vector<double> writeVectorData(meshSize * dim);
      std::vector<double> justInTimeWriteVectorData(meshSize * dim);

      if (time == 1) {
        // Linear filling for writeData1
        std::generate(writeData.begin(), writeData.end(), [n = 0]() mutable {
          return 0.8 * n++; // Simple increasing pattern, starting from 0
        });
        // vectorData1 with linearly increasing values
        for (int i = 0; i < meshSize; ++i) {
          writeVectorData[dim * i]     = 12.5 * i;  // x component
          writeVectorData[dim * i + 1] = 12.2 * i;  // y component
          writeVectorData[dim * i + 2] = -12.2 * i; // z component
        }
      } else if (time == 2) {
        // Quadratic filling for writeData2
        std::generate(writeData.begin(), writeData.end(), [n = 0]() mutable {
          return -28.0 + 0.1 * (n * n++); // Quadratic pattern, n^2 scaled by 0.1
        });
        // vectorData2 with quadratic and linear patterns
        for (int i = 0; i < meshSize; ++i) {
          writeVectorData[dim * i]     = 0.28 * i * i;          // x component quadratically increasing
          writeVectorData[dim * i + 1] = 100 - 0.5 * i;         // y component linearly decreasing
          writeVectorData[dim * i + 2] = 100 - 0.5 * i * i * i; // z component linearly decreasing
        }
      } else {
        PRECICE_ASSERT(false);
      }
      // First, we check the separate polynomial PUM (scalar and vector)
      couplingInterface.writeData(meshName, dataName, ids, writeData);
      couplingInterface.writeAndMapData(otherMeshName, dataName, positions, writeData);
      couplingInterface.writeData(meshName, vectorDataName, ids, writeVectorData);
      couplingInterface.writeAndMapData(otherMeshName, vectorDataName, positions, writeVectorData);

      // Second, we check the no polynomial PUM (scalar and vector)
      couplingInterface.writeAndMapData(otherMeshBName, dataName, positions, writeData);
      couplingInterface.writeAndMapData(otherMeshBName, vectorDataName, positions, writeVectorData);

      // solve time step
      // write data (not necessary here)
      couplingInterface.advance(dt);
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName       = "MeshTwo";
    auto meshBName      = "MeshTwoB";
    auto meshCName      = "MeshTwoC";
    auto meshDName      = "MeshTwoD";
    auto dataName       = "Temperature";
    auto vectorDataName = "Velocity";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName));

    std::vector<double> positions;
    int                 size = 250;
    positions.reserve(size * dim);
    for (int i = 0; i < 10; ++i) {
      for (int j = 0; j < 5; ++j) {
        for (int k = 0; k < 5; ++k) {
          positions.emplace_back(i * 0.1);
          positions.emplace_back(j * 0.2);
          positions.emplace_back(k * 0.2);
        }
      }
    }
    std::vector<int> idsA(size, -1);
    std::vector<int> idsB(size, -1);
    std::vector<int> idsC(size, -1);
    std::vector<int> idsD(size, -1);
    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, idsA);
    couplingInterface.setMeshVertices(meshBName, positions, idsB);
    couplingInterface.setMeshVertices(meshCName, positions, idsC);
    couplingInterface.setMeshVertices(meshDName, positions, idsD);
    // Some dummy readData
    std::vector<double> readData(size);
    std::vector<double> justInTimeReadData(size);
    std::vector<double> readVectorData(size * dim);
    std::vector<double> justInTimeReadVectorData(size * dim);

    // Initialize
    couplingInterface.initialize();
    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // Meshes A and C have the same mapping
      // just-in-time mapping
      couplingInterface.readData(meshName, dataName, idsA, dt, justInTimeReadData);
      couplingInterface.readData(meshName, vectorDataName, idsA, dt, justInTimeReadVectorData);
      // conventional
      couplingInterface.readData(meshCName, dataName, idsC, dt, readData);
      couplingInterface.readData(meshCName, vectorDataName, idsC, dt, readVectorData);

      BOOST_TEST(justInTimeReadData == readData, boost::test_tools::per_element());
      BOOST_TEST(justInTimeReadVectorData == readVectorData, boost::test_tools::per_element());

      // Meshes B and D have the same mapping
      // just-in-time mapping
      couplingInterface.readData(meshBName, dataName, idsB, dt, justInTimeReadData);
      couplingInterface.readData(meshBName, vectorDataName, idsB, dt, justInTimeReadVectorData);
      // conventional
      couplingInterface.readData(meshDName, dataName, idsD, dt, readData);
      couplingInterface.readData(meshDName, vectorDataName, idsD, dt, readVectorData);

      BOOST_TEST(justInTimeReadData == readData, boost::test_tools::per_element());
      BOOST_TEST(justInTimeReadVectorData == readVectorData, boost::test_tools::per_element());

      // read data
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
