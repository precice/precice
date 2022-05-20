#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(testParallelSquare1To2)
{
  using precice::VertexID;
  using precice::testing::equals;

  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(2_ranks));

  // Implement your test here.
  BOOST_TEST(true);
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;

  if (context.isNamed("SolverOne")) {
    auto meshID = interface.getMeshID("MeshOne");
    auto dataID = interface.getDataID("DataOne", meshID);

    std::vector<double> coords;

    // Create a square with top left corner (rank 0) or bottom right. Diagonal "y = x" is shared.
    coords = {0.0, 0.0,
              1.0, 0.0,
              1.0, 1.0,
              0.0, 1.0};
    vertexIDs.resize(coords.size() / 2);
    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());

    // Square ABCD in counter-clockwise order. A is the origin, B on the right
    auto AB = interface.setMeshEdge(meshID, vertexIDs[0], vertexIDs[1]);
    auto BC = interface.setMeshEdge(meshID, vertexIDs[1], vertexIDs[2]);
    auto CD = interface.setMeshEdge(meshID, vertexIDs[2], vertexIDs[3]);
    auto DA = interface.setMeshEdge(meshID, vertexIDs[3], vertexIDs[0]);
    auto CA = interface.setMeshEdge(meshID, vertexIDs[2], vertexIDs[0]);

    interface.setMeshTriangle(meshID, AB, BC, CA);
    interface.setMeshTriangle(meshID, CA, CD, DA);

    dt = interface.initialize();

    // Run a step and write data with f(x) = x+2*y
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    values = {0.0,
              1.0,
              3.0,
              2.0};

    interface.writeBlockScalarData(dataID, 4, vertexIDs.data(), values.data());

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();
  } else { // SolverTwo
    auto meshID = interface.getMeshID("MeshTwo");
    auto dataID = interface.getDataID("DataOne", meshID);

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
    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());

    dt = interface.initialize();

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

    interface.readBlockScalarData(dataID, expected.size(), vertexIDs.data(), readData.data());
    BOOST_CHECK(equals(expected, readData));
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
