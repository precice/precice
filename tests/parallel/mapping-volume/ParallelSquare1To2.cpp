#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(ParallelSquare1To2)
{
  using precice::VertexID;
  using precice::testing::equals;

  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(2_ranks));

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;

  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne";

    std::vector<double> coords;

    // Create a square with top left corner (rank 0) or bottom right. Diagonal "y = x" is shared.
    coords = {0.0, 0.0,
              1.0, 0.0,
              1.0, 1.0,
              0.0, 1.0};
    vertexIDs.resize(coords.size() / 2);
    interface.setMeshVertices(meshName, vertexIDs.size(), coords.data(), vertexIDs.data());

    // Square ABCD in counter-clockwise order. A is the origin, B on the right
    interface.setMeshTriangle(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2]);
    interface.setMeshTriangle(meshName, vertexIDs[0], vertexIDs[2], vertexIDs[3]);

    interface.initialize();
    dt = interface.getMaxTimeStepSize();

    // Run a step and write data with f(x) = x+2*y
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    values = {0.0,
              1.0,
              3.0,
              2.0};

    interface.writeBlockScalarData(meshName, dataName, 4, vertexIDs.data(), values.data());

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();
  } else { // SolverTwo
    auto meshName = "MeshTwo";
    auto dataName = "DataOne";

    std::vector<double> coords;
    if (context.rank == 0) {
      coords = {1. / 6, 1. / 2,
                1. / 2, 1. / 6};
    } else {
      coords = {
          5. / 6, 1. / 2,
          1. / 2, 5. / 6};
    }

    vertexIDs.resize(coords.size() / 2);
    interface.setMeshVertices(meshName, vertexIDs.size(), coords.data(), vertexIDs.data());

    interface.initialize();
    dt = interface.getMaxTimeStepSize();

    // Run a step and read data expected to be f(x) = x+2*y
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(2);
    Eigen::VectorXd readData(2);
    if (context.rank == 0) {
      expected << 7. / 6, 5. / 6;
    } else {
      expected << 11. / 6, 13. / 6;
    }

    dt = interface.getMaxTimeStepSize();
    interface.readBlockScalarData(meshName, dataName, expected.size(), vertexIDs.data(), dt, readData.data());
    BOOST_CHECK(equals(expected, readData));
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
