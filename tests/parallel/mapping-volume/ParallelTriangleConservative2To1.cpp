#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(ParallelTriangleConservative2To1)
{
  using precice::VertexID;
  using precice::testing::equals;

  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(1_rank));
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;

  if (context.isNamed("SolverOne")) {
    auto meshID = "MeshOne";
    auto dataID = "DataOne"; //  meshID

    std::vector<double> coords;

    // Each rank sends one "force" on one point
    if (context.rank == 0) {
      coords = {0.3, 0.5};
    } else {
      coords = {0.7, 0.2};
    }
    vertexIDs.resize(coords.size() / 2);
    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());

    dt = interface.initialize();

    // Run a step and write forces
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    values = {1.0};

    interface.writeBlockScalarData(dataID, 1, vertexIDs.data(), values.data());

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();
  } else { // SolverTwo
    auto meshID = "MeshTwo";
    auto dataID = "DataOne"; //  meshID

    std::vector<double> coords = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0}; // Lower-left triangle making half the unit square

    vertexIDs.resize(coords.size() / 2);
    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());
    interface.setMeshTriangle(meshID, vertexIDs[0], vertexIDs[1], vertexIDs[2]);

    dt = interface.initialize();

    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(3);
    Eigen::VectorXd readData(3);
    // Sum of two forces: one is {0.2, 0.3, 0.5} and the second is {0.1, 0.7, 0.2}
    // These are proportional to barycentric coordinates.
    expected << 0.3, 1.0, 0.7;

    interface.readBlockScalarData(dataID, expected.size(), vertexIDs.data(), readData.data());
    BOOST_CHECK(equals(expected, readData));
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
