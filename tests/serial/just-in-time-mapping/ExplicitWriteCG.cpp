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
// cg-conservative-write

BOOST_AUTO_TEST_CASE(ExplicitWriteCG)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 3;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};

  if (context.isNamed("SolverOne")) {
    auto otherMeshName  = "MeshTwo";
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
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);

    couplingInterface.initialize();

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // write data
      std::vector<double> writeVectorData(meshSize * dim);

      if (time == 1) {
        // vectorData1 with linearly increasing values
        for (int i = 0; i < meshSize; ++i) {
          writeVectorData[dim * i]     = 0.1 * i; // x component
          writeVectorData[dim * i + 1] = 0.2 * i; // y component
          writeVectorData[dim * i + 2] = 0.4 * i; // z component
        }
      } else if (time == 2) {
        // Quadratic filling for writeData2
        // vectorData2 with quadratic and linear patterns
        for (int i = 0; i < meshSize; ++i) {
          writeVectorData[dim * i]     = 0; // x component quadratically increasing
          writeVectorData[dim * i + 1] = 0; // y component linearly decreasing
          writeVectorData[dim * i + 2] = 0; // z component linearly decreasing
        }
      } else {
        BOOST_TEST(false);
      }
      couplingInterface.writeAndMapData(otherMeshName, vectorDataName, positions, writeVectorData);

      // solve time step
      // write data (not necessary here)
      couplingInterface.advance(dt);
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName       = "MeshTwo";
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

    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, idsA);
    // Some dummy readData
    std::vector<double> justInTimeReadVectorData(size * dim);

    // Initialize
    couplingInterface.initialize();
    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;
      // Meshes A and C have the same mapping
      // just-in-time mapping
      couplingInterface.readData(meshName, vectorDataName, idsA, dt, justInTimeReadVectorData);

      if (time == 1) {
        BOOST_TEST(justInTimeReadVectorData[0] == 0);
        BOOST_TEST(justInTimeReadVectorData[1] == 0);
        BOOST_TEST(justInTimeReadVectorData[2] == 0);
        double sample = 15.44379927617392;
        BOOST_TEST(justInTimeReadVectorData[3] == sample);
        BOOST_TEST(justInTimeReadVectorData[4] == 2 * sample);
        BOOST_TEST(justInTimeReadVectorData[5] == 4 * sample);
      } else {
        Eigen::VectorXd expected(size * dim);
        expected.setZero();
        BOOST_TEST(justInTimeReadVectorData == expected, boost::test_tools::per_element());
      }
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
