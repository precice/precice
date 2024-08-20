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

  // For velocity
  std::vector<double> velocityData1({27, 12, 8});
  std::vector<double> expectedVelocityData1({27, 20, 0, 0});

  // For vorticity
  std::vector<double> vorticityData1({8000, 4, 19, 5.5, 8, 9, 0.5, 7, 9});
  std::vector<double> expectedVorticityData1({8000, 4, 19, 6, 15, 18, 0, 0, 0, 0, 0, 0});

  // for porosity
  std::vector<double> readPorosityData1(3, -1);
  std::vector<double> porosityData1({100, 500, 700, 900});
  std::vector<double> expectedPorosityData1({0, 0, 0}); // no initial data

  // for speed of light
  std::vector<double> readSpeedOfLightData1(3 * dim, -1);
  std::vector<double> SpeedOfLightData1({0.3, -4, 7, 28, 0, 4, 0.03, 0.04, 0.05, 1, 4, 7});
  std::vector<double> expectedSpeedOfLightData1({0, 0, 0, 0, 0, 0, 0, 0, 0}); // no initial data

  // maps to 2, 2, 3, 0 of the output mesh
  std::vector<double> tmpPositions2 = {0.5, 0.45, 0.6, 0.5, 0.55, 0.6, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

  std::vector<double> velocityData2({4, -7, 28, 3500});
  std::vector<double> vorticityData2({1000, 2000, 3000, 4, 5, 8, 0.1, 0.2, 0.3, 15., 3, 8});
  std::vector<double> expectedVelocityData2({3500, 0, -3, 28});
  std::vector<double> expectedVorticityData2({15., 3, 8, 0, 0, 0, 1004, 2005, 3008, 0.1, 0.2, 0.3});
  std::vector<double> readPorosityData2(4, -1);
  std::vector<double> expectedPorosityData2({700, 700, 900, 100});
  std::vector<double> readSpeedOfLightData2(4 * dim, -1);
  std::vector<double> expectedSpeedOfLightData2({0.03, 0.04, 0.05, 0.03, 0.04, 0.05, 1, 4, 7, 0.3, -4, 7});

  auto velocityData     = "Velocities";
  auto vorticityData    = "Vorticity";
  auto porosityData     = "Porosity";
  auto speedOfLightData = "SpeefOfLight";
  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 3);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);

    couplingInterface.initialize();

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // read data
      if (time == 1) {
        for (std::size_t i = 0; i < readPorosityData1.size(); ++i) {
          std::vector<double> solverTwoCoord(dim);
          for (int d = 0; d < dim; ++d) {
            solverTwoCoord[d] = tmpPositions1[i * dim + d];
          }
          couplingInterface.mapAndreadData(otherMeshName, porosityData, solverTwoCoord, dt, {&readPorosityData1[i], 1});
          couplingInterface.mapAndreadData(otherMeshName, speedOfLightData, solverTwoCoord, dt, {&readSpeedOfLightData1[i * dim], dim});
        }

        BOOST_TEST(readPorosityData1 == expectedPorosityData1);
        BOOST_TEST(readSpeedOfLightData1 == expectedSpeedOfLightData1);

      } else if (time == 2) {
        couplingInterface.mapAndreadData(otherMeshName, porosityData, tmpPositions2, dt, readPorosityData2);
        BOOST_TEST(readPorosityData2 == expectedPorosityData2);

        couplingInterface.mapAndreadData(otherMeshName, speedOfLightData, tmpPositions2, dt, readSpeedOfLightData2);
        BOOST_TEST(readSpeedOfLightData2 == expectedSpeedOfLightData2);
      } else {
        PRECICE_ASSERT(false);
      }
      // solve time step
      // write data:
      if (time == 1) {
        for (std::size_t i = 0; i < velocityData1.size(); ++i) {
          std::vector<double> solverTwoCoord(dim);
          for (int d = 0; d < dim; ++d) {
            solverTwoCoord[d] = tmpPositions1[i * dim + d];
          }
          couplingInterface.mapAndwriteData(otherMeshName, velocityData, solverTwoCoord, {&velocityData1[i], 1});
          couplingInterface.mapAndwriteData(otherMeshName, vorticityData, solverTwoCoord, {&vorticityData1[i * dim], dim});
        }
      } else if (time == 2) {
        // The second time, we pass all data at once
        couplingInterface.mapAndwriteData(otherMeshName, velocityData, tmpPositions2, velocityData2);
        couplingInterface.mapAndwriteData(otherMeshName, vorticityData, tmpPositions2, vorticityData2);
      } else {
        PRECICE_ASSERT(false);
      }
      couplingInterface.advance(dt);
    }

    // @todo: this call is allowed. Which data should it return? for dt = 0, it should correspond to the previous dt = 1, right?
    // the last data written by the second participant is essentially never received on the other side
    //  couplingInterface.mapAndreadData(otherMeshName, porosityData, tmpPositions2, 0, readPorosityData2);

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName = "MeshTwo";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName) == 3);

    std::vector<double> positions = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};
    std::vector<int>    ids(4, -1);

    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // Some dummy readData
    std::vector<double> readData(4, -1);
    std::vector<double> readData2(4 * dim, -1);

    double time = 0;
    // Initialize
    couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // read data
      couplingInterface.readData(meshName, velocityData, ids, dt, readData);
      couplingInterface.readData(meshName, vorticityData, ids, dt, readData2);
      if (time == 1) {
        BOOST_TEST(std::accumulate(readData.begin(), readData.end(), 0.) == std::accumulate(velocityData1.begin(), velocityData1.end(), 0.));
        BOOST_TEST(readData == expectedVelocityData1);
        BOOST_TEST(expectedVorticityData1 == readData2);
      } else if (time == 2) {
        BOOST_TEST(std::accumulate(readData.begin(), readData.end(), 0.) == std::accumulate(velocityData2.begin(), velocityData2.end(), 0.));
        BOOST_TEST(readData == expectedVelocityData2);
        BOOST_TEST(expectedVorticityData2 == readData2);
      } else {
        PRECICE_ASSERT(false);
      }
      // solve time step
      //write data
      couplingInterface.writeData(meshName, porosityData, ids, porosityData1);
      couplingInterface.writeData(meshName, speedOfLightData, ids, SpeedOfLightData1);

      couplingInterface.advance(dt);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
