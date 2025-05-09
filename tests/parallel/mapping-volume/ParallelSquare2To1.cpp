#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
PRECICE_TEST_SETUP("SolverOne"_on(2_ranks), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ParallelSquare2To1)
{
  PRECICE_TEST();

  using precice::VertexID;
  using precice::testing::equals;

  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;

  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne";

    std::vector<double> coords;

    // Create a square with top left corner (rank 0) or bottom right. Diagonal "y = x" is shared.
    if (context.rank == 0) {
      coords = {0.0, 0.0,
                1.0, 1.0,
                0.0, 1.0};
    } else {
      coords = {0.0, 0.0,
                1.0, 1.0,
                1.0, 0.0};
    }

    vertexIDs.resize(coords.size() / 2);
    interface.setMeshVertices(meshName, coords, vertexIDs);
    interface.setMeshTriangle(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2]);

    interface.initialize();
    dt = interface.getMaxTimeStepSize();

    // Run a step and write data with f(x) = x+2*y
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    if (context.rank == 0) {
      values = {0.0,
                3.0,
                2.0};
    } else {
      values = {0.0,
                3.0,
                1.0};
    }

    interface.writeData(meshName, dataName, vertexIDs, values);

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();

  } else {
    auto meshName = "MeshTwo";
    auto dataName = "DataOne";

    std::vector<double> coords;

    /* Each ranks reads points from both input meshes
    Rank 0: (1/6, 1/2) and (1/2, 1/6)
    Rank 1: (5/6, 1/2) and (1/2, 5/6)

    */
    coords = {1. / 6, 1. / 2,
              1. / 2, 1. / 6,
              5. / 6, 1. / 2,
              1. / 2, 5. / 6};

    vertexIDs.resize(coords.size() / 2);
    interface.setMeshVertices(meshName, coords, vertexIDs);

    interface.initialize();
    dt = interface.getMaxTimeStepSize();

    // Run a step and read data expected to be f(x) = x+2*y
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(4);
    Eigen::VectorXd readData(4);
    expected << 7. / 6, 5. / 6, 11. / 6, 13. / 6;

    dt = interface.getMaxTimeStepSize();
    interface.readData(meshName, dataName, vertexIDs, dt, readData);
    BOOST_CHECK(equals(expected, readData));
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
