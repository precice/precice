#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(ParallelCube3To1)
{
  PRECICE_TEST("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank));
  // Split the unit cube in 6 tetrahedra (2 per rank) and set up a consistent mapping

  using precice::VertexID;
  using precice::testing::equals;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;
  double                a = 1, b = 2, c = 5, d = 1;

  if (context.isNamed("SolverOne")) {
    auto meshID = "MeshOne";
    auto dataID = "DataOne"; //  meshID

    std::vector<double> coords;

    switch (context.rank) {
    case 0:
      coords = {0, 0, 0,
                0, 0, 1,
                0, 1, 1,
                1, 1, 1,
                0, 1, 0};
      break;
    case 1:
      coords = {0, 0, 0,
                0, 0, 1,
                1, 0, 1,
                1, 1, 1,
                1, 0, 0};
      break;
    case 2:
      coords = {0, 0, 0,
                0, 1, 0,
                1, 1, 1,
                1, 1, 0,
                1, 0, 0};
      break;
    }

    vertexIDs.resize(coords.size() / 3);
    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());
    switch (context.rank) {
    case 0:
      interface.setMeshTetrahedron(meshID, vertexIDs[0], vertexIDs[1], vertexIDs[2], vertexIDs[3]);
      interface.setMeshTetrahedron(meshID, vertexIDs[0], vertexIDs[4], vertexIDs[2], vertexIDs[3]);
      break;
    case 1:
      interface.setMeshTetrahedron(meshID, vertexIDs[0], vertexIDs[1], vertexIDs[2], vertexIDs[3]);
      interface.setMeshTetrahedron(meshID, vertexIDs[0], vertexIDs[4], vertexIDs[2], vertexIDs[3]);
      break;
    case 2:
      interface.setMeshTetrahedron(meshID, vertexIDs[0], vertexIDs[1], vertexIDs[2], vertexIDs[3]);
      interface.setMeshTetrahedron(meshID, vertexIDs[0], vertexIDs[4], vertexIDs[2], vertexIDs[3]);
      break;
    }

    dt = interface.initialize();

    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    // Sample the function ax + by + cz + d
    for (std::size_t i = 0; i < vertexIDs.size(); ++i) {
      values.push_back(d + a * coords[3 * i] + b * coords[3 * i + 1] + c * coords[3 * i + 2]);
    }

    interface.writeBlockScalarData(meshID, dataID, vertexIDs.size(), vertexIDs.data(), values.data());

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();

  } else {
    auto meshID = "MeshTwo";
    auto dataID = "DataOne"; //  meshID

    // For completion, we sample points in each direction.
    std::vector<double> coords;
    std::vector<double> values;
    const int           SAMPLING = 5;
    for (int i = 0; i <= SAMPLING; ++i) {
      for (int j = 0; j <= SAMPLING; ++j) {
        for (int k = 0; k <= SAMPLING; ++k) {
          double x = double(i) / SAMPLING;
          double y = double(j) / SAMPLING;
          double z = double(k) / SAMPLING;
          coords.push_back(x);
          coords.push_back(y);
          coords.push_back(z);
          values.push_back(d + a * x + b * y + c * z);
        }
      }
    }

    vertexIDs.resize(coords.size() / 3);
    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());

    dt = interface.initialize();

    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    //Check expected VS read
    Eigen::VectorXd expected(values.size());
    for (int i = 0; i < expected.size(); ++i) {
      expected(i) = values[i];
    }
    Eigen::VectorXd readData(values.size());

    interface.readBlockScalarData(meshID, dataID, expected.size(), vertexIDs.data(), readData.data());
    //BOOST_CHECK(equals(expected, readData, 1e-3));
    for (int i = 0; i < expected.size(); ++i) {
      BOOST_CHECK(equals(expected(i), readData[i]));
    }
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
