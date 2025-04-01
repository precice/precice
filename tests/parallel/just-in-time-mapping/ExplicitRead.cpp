#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(JustInTimeMapping)
PRECICE_TEST_SETUP("SolverOne"_on(2_ranks), "SolverTwo"_on(1_rank))

// Test case for a just-in-time mapping on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. SolverTwo defines the mesh (as usual) whereas SolverOne maps
// just-in-time from this mesh.
// parallel-access-nearest-neighbor-consistent-read
BOOST_AUTO_TEST_CASE(ExplicitRead)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), context.rank, context.size);

  constexpr int dim = 2;

  if (context.isNamed("SolverOne")) {

    std::array<double, dim * 2> boundingBox;
    std::vector<double>         tmpPositions;
    std::vector<double>         expectedData1;
    std::vector<double>         expectedData2;

    if (context.isPrimary()) {
      boundingBox   = {0.0, 1.0, 0.0, 1.0};
      tmpPositions  = {0.1, 0.1, 0.1, 0.5, 1.0, 0.0, 1.0, 1.0};
      expectedData1 = {1, 1, 2, 4};
      expectedData2 = {-10, -10, -11, -13};
    }
    if (!context.isPrimary()) {
      boundingBox   = {1.0, 2.0, 0.0, 1.0};
      tmpPositions  = {1.2, 0.0, 2.0, 1.0, 1.0, 0.5};
      expectedData1 = {2, 5, 3};
      expectedData2 = {-11, -15, -12};
    }
    auto otherMeshName = "MeshTwo";
    auto dataName      = "Velocities";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 2);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);

    couplingInterface.initialize();

    // Although we don't use the direct mesh access here, we check that the initialization works
    // the repartitioned rank-local mesh is here filered according to the access region
    // the local mesh contains for both ranks more vertices
    {
      auto meshSize = couplingInterface.getMeshVertexSize(otherMeshName);
      BOOST_TEST(meshSize == (context.isPrimary() ? 2 : 3));
      std::vector<double> receivedCoordinates(meshSize * dim);
      std::vector<int>    receivedIDs(meshSize);
      couplingInterface.getMeshVertexIDsAndCoordinates(otherMeshName, receivedIDs, receivedCoordinates);
      std::vector<double> expectedCoordinates = context.isPrimary() ? std::vector<double>({0.0, 0.0, 0.9, 0.8}) : std::vector<double>({1.2, 0.1, 1.2, 1.0, 2.0, 1.0});
      BOOST_TEST(receivedCoordinates == expectedCoordinates, boost::test_tools::per_element());
    }

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

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
          BOOST_TEST(false);
        }
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

    std::vector<double> positions = {0.0, 0.0, 1.2, 0.1, 0.9, 0.8, 1.2, 1.0, 2.0, 1.0};
    std::vector<int>    ids(5, -1);

    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // Some dummy readData
    std::vector<double> writeData1({1, 2, 3, 4, 5});
    std::vector<double> writeData2({-10, -11, -12, -13, -15});

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
