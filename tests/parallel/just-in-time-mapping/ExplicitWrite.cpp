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
// parallel-access-nearest-neighbor-conservative-write
BOOST_AUTO_TEST_CASE(ExplicitWrite)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), context.rank, context.size);

  constexpr int dim = 2;

  if (context.isNamed("SolverOne")) {

    std::array<double, dim * 2> boundingBox;
    std::vector<double>         tmpPositions;

    if (context.isPrimary()) {
      boundingBox  = {0.0, 1.0, 0.0, 1.0};
      tmpPositions = {0.1, 0.1, 0.1, 0.5, 1.0, 0.0, 1.0, 1.0};
    }
    if (!context.isPrimary()) {
      boundingBox  = {1.0, 2.0, 0.0, 1.0};
      tmpPositions = {1.2, 0.0, 2.0, 1.0, 1.0, 0.5};
    }
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

      for (std::size_t i = 0; i < tmpPositions.size() / dim; ++i) {
        std::vector<double> solverTwoCoord(dim);
        double              value = time * 10 + i - 7 * context.rank;
        for (int d = 0; d < dim; ++d) {
          solverTwoCoord[d] = tmpPositions[i * dim + d];
        }
        couplingInterface.mapAndWriteData(otherMeshName, dataName, solverTwoCoord, {&value, 1});
      }
      // solve time step
      // write data (not necessary here)
      couplingInterface.advance(dt);
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));

    std::vector<double> expectedData1 = {21, 15, 5, 13, 4};
    std::vector<double> expectedData2 = {41, 35, 15, 23, 14};

    // Query IDs
    auto meshName = "MeshTwo";
    auto dataName = "Velocities";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName));

    std::vector<double> positions = {0.0, 0.0, 1.2, 0.1, 0.9, 0.8, 1.2, 1.0, 2.0, 1.0};
    std::vector<int>    ids(5, -1);

    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // Some dummy readData
    std::vector<double> readData(5, -1);
    // Initialize
    couplingInterface.initialize();
    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;
      // read data
      couplingInterface.readData(meshName, dataName, ids, dt, readData);

      if (time == 1) {
        BOOST_TEST(expectedData1 == readData, boost::test_tools::per_element());
      } else if (time == 2) {
        BOOST_TEST(expectedData2 == readData, boost::test_tools::per_element());
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
