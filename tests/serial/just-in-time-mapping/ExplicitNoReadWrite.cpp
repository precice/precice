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
// by another participant.
// Here, we use a just-in-time mapping, but don't call the read and write functions
// in the first time step (which should be ok)
BOOST_AUTO_TEST_CASE(ExplicitNoReadWrite)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 3;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};

  auto velocityData = "Velocities";
  auto porosityData = "Porosity";
  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 3);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);

    couplingInterface.initialize();

    std::vector<double> positions;
    int                 size = 20;
    positions.reserve(size * dim);
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
          positions.emplace_back(i * 0.2);
          positions.emplace_back(j * 0.5);
          positions.emplace_back(k * 0.5);
        }
      }
    }
    std::vector<double> readData(size, -1);
    std::vector<double> ref(size, 5);

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // read data
      if (time == 1) {
        // Do nothing
      } else if (time == 2) {
        couplingInterface.mapAndReadData(otherMeshName, porosityData, positions, dt, readData);
        // BOOST_TEST(std::all_of(readData.begin(), readData.end(), [](auto r) { return r == 5.0; }));
        BOOST_TEST(readData == ref, boost::test_tools::per_element());
      } else {
        BOOST_TEST(false);
      }
      // solve time step
      // write data:
      if (time == 1) {
        // Do nothing
      } else if (time == 2) {
        // The second time, we pass all data at once
        couplingInterface.writeAndMapData(otherMeshName, velocityData, positions, ref);
      } else {
        BOOST_TEST(false);
      }
      couplingInterface.advance(dt);
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName = "MeshTwo";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName) == 3);

    std::vector<double> positions;
    int                 size = 30;
    positions.reserve(size * dim);
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 2; ++k) {
          positions.emplace_back(i * 0.2);
          positions.emplace_back(j * 0.3);
          positions.emplace_back(k * 0.5);
        }
      }
    }

    std::vector<int> ids(30, -1);
    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // Some dummy readData
    std::vector<double> readData(size, -1);

    double time = 0;
    // Initialize
    couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // read data
      couplingInterface.readData(meshName, velocityData, ids, dt, readData);
      if (time == 1) {
        BOOST_TEST(std::accumulate(readData.begin(), readData.end(), 0.) == 0);
      } else if (time == 2) {
        // Just do a cumulative check here (20 x 5 = vertices x val written at t = 2)
        // It's mostly to check that the caches are allocated correctly
        BOOST_TEST(std::accumulate(readData.begin(), readData.end(), 0.) == (20 * 5));
      } else {
        BOOST_TEST(false);
      }
      // solve time step
      // write data
      std::vector<double> writeData(size, 5);
      couplingInterface.writeData(meshName, porosityData, ids, writeData);

      couplingInterface.advance(dt);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
