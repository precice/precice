#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
PRECICE_TEST_SETUP("SolverOne"_on(3_ranks), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ParallelCube3To1)
{
  PRECICE_TEST();
  // Split the unit cube in 6 tetrahedra (2 per rank) and set up a consistent mapping

  using precice::VertexID;
  using precice::testing::equals;

  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;
  double                a = 1, b = 2, c = 5, d = 1;

  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne";

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
    interface.setMeshVertices(meshName, coords, vertexIDs);
    switch (context.rank) {
    case 0:
      interface.setMeshTetrahedron(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2], vertexIDs[3]);
      interface.setMeshTetrahedron(meshName, vertexIDs[0], vertexIDs[4], vertexIDs[2], vertexIDs[3]);
      break;
    case 1:
      interface.setMeshTetrahedron(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2], vertexIDs[3]);
      interface.setMeshTetrahedron(meshName, vertexIDs[0], vertexIDs[4], vertexIDs[2], vertexIDs[3]);
      break;
    case 2:
      interface.setMeshTetrahedron(meshName, vertexIDs[0], vertexIDs[1], vertexIDs[2], vertexIDs[3]);
      interface.setMeshTetrahedron(meshName, vertexIDs[0], vertexIDs[4], vertexIDs[2], vertexIDs[3]);
      break;
    }

    interface.initialize();
    dt = interface.getMaxTimeStepSize();

    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    // Sample the function ax + by + cz + d
    for (std::size_t i = 0; i < vertexIDs.size(); ++i) {
      values.push_back(d + a * coords[3 * i] + b * coords[3 * i + 1] + c * coords[3 * i + 2]);
    }

    interface.writeData(meshName, dataName, vertexIDs, values);

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();

  } else {
    auto meshName = "MeshTwo";
    auto dataName = "DataOne";

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
    interface.setMeshVertices(meshName, coords, vertexIDs);

    interface.initialize();

    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");
    dt = interface.getMaxTimeStepSize();
    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(values.size());
    for (int i = 0; i < expected.size(); ++i) {
      expected(i) = values[i];
    }
    Eigen::VectorXd readData(values.size());
    dt = interface.getMaxTimeStepSize();
    interface.readData(meshName, dataName, vertexIDs, dt, readData);
    // BOOST_CHECK(equals(expected, readData, 1e-3));
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
