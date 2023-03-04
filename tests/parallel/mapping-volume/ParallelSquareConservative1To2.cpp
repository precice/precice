#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(ParallelSquareConservative1To2)
{
  using precice::VertexID;
  using precice::testing::equals;

  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(2_ranks));
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;

  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne"; //  meshName

    std::vector<double> coords;

    // Apply forces on two points
    coords = {0.3, 0.5,
              0.9, 0.2};
    vertexIDs.resize(coords.size() / 2);
    interface.setMeshVertices(meshName, vertexIDs.size(), coords.data(), vertexIDs.data());

    dt = interface.initialize();

    // Run a step and write forces
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    values = {1.0,
              1.0};

    interface.writeBlockScalarData(meshName, dataName, 2, vertexIDs.data(), values.data());

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();
  } else { // SolverTwo
    auto meshName = "MeshTwo";
    auto dataName = "DataOne"; //  meshName

    std::vector<double> coords;
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
    interface.setMeshVertices(meshName, vertexIDs.size(), coords.data(), vertexIDs.data());
    interface.setMeshTriangle(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2]);

    dt = interface.initialize();

    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(3);
    Eigen::VectorXd readData(3);
    if (context.rank == 0) {
      // Force applied on (0.3, 0.5) to spread with respect to barycentric coordinates
      expected << 0.5, 0.3, 0.2;
    } else {
      // FOrce applied on (0.9, 0.2)
      expected << 0.1, 0.2, 0.7;
    }

    interface.readBlockScalarData(meshName, dataName, expected.size(), vertexIDs.data(), readData.data());
    BOOST_CHECK(equals(expected, readData));
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
