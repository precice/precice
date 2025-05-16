#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <iostream>
#include <precice/precice.hpp>
#include <vector>
#include "precice/impl/ParticipantImpl.hpp"
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(3_ranks))
BOOST_AUTO_TEST_CASE(ParallelCube1To3)
{
  PRECICE_TEST();
  using precice::VertexID;
  using precice::testing::equals;

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;
  double                a = 1, b = 2, c = 5, d = 1;

  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne";

    std::vector<double> coords;

    // Unit cube made of 6 tetra, with 8 points
    coords = {0, 0, 0,
              0, 0, 1,
              0, 1, 0,
              0, 1, 1,
              1, 0, 0,
              1, 0, 1,
              1, 1, 0,
              1, 1, 1};

    vertexIDs.resize(coords.size() / 3);
    participant.setMeshVertices(meshName, coords, vertexIDs);

    VertexID v000 = vertexIDs[0];
    VertexID v001 = vertexIDs[1];
    VertexID v010 = vertexIDs[2];
    VertexID v011 = vertexIDs[3];
    VertexID v100 = vertexIDs[4];
    VertexID v101 = vertexIDs[5];
    VertexID v110 = vertexIDs[6];
    VertexID v111 = vertexIDs[7];

    participant.setMeshTetrahedron(meshName, v000, v001, v011, v111);
    participant.setMeshTetrahedron(meshName, v000, v010, v011, v111);
    participant.setMeshTetrahedron(meshName, v000, v001, v101, v111);
    participant.setMeshTetrahedron(meshName, v000, v100, v101, v111);
    participant.setMeshTetrahedron(meshName, v000, v010, v110, v111);
    participant.setMeshTetrahedron(meshName, v000, v100, v110, v111);

    auto &mesh = precice::testing::WhiteboxAccessor::impl(participant).mesh("MeshOne");
    BOOST_REQUIRE(mesh.nVertices() == 8);
    BOOST_REQUIRE(mesh.tetrahedra().size() == 6);

    participant.initialize();
    dt = participant.getMaxTimeStepSize();

    // Run a step and write data with f(x) = ax + by + cz + d
    BOOST_TEST(participant.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    values = {d,
              c + d,
              b + d,
              b + c + d,
              a + d,
              a + c + d,
              a + b + d,
              c + a + b + d};

    participant.writeData(meshName, dataName, vertexIDs, values);

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Sending participant must advance only once.");
    participant.finalize();
  } else { // SolverTwo
    auto meshName = "MeshTwo";
    auto dataName = "DataOne";

    std::vector<double> coords;

    const int SAMPLING = 5;

    switch (context.rank) {
    case 0:
      // Sample the sub-region "x,y and z > 0.5"
      for (int i = 0; i <= SAMPLING; ++i) {
        for (int j = 0; j <= SAMPLING; ++j) {
          for (int k = 0; k <= SAMPLING; ++k) {
            double x = 0.5 + 0.5 * double(i) / SAMPLING;
            double y = double(j) / SAMPLING;
            double z = double(k) / SAMPLING;
            coords.push_back(x);
            coords.push_back(y);
            coords.push_back(z);
          }
        }
      }
      break;
    case 1:
      // Sample the sub-region "x > 0.5, y < 0.5, all z"
      for (int i = 0; i <= SAMPLING; ++i) {
        for (int j = 0; j <= SAMPLING; ++j) {
          for (int k = 0; k <= SAMPLING; ++k) {
            double x = 0.5 + 0.5 * double(i) / SAMPLING;
            double y = 0.5 * double(j) / SAMPLING;
            double z = double(k) / SAMPLING;
            coords.push_back(x);
            coords.push_back(y);
            coords.push_back(z);
          }
        }
      }
      break;
    case 2:
      // Sample the center:  "x, y and z between 0.4 and 0.6"
      for (int i = 0; i <= SAMPLING; ++i) {
        for (int j = 0; j <= SAMPLING; ++j) {
          for (int k = 0; k <= SAMPLING; ++k) {
            double x = 0.4 + 0.2 * double(i) / SAMPLING;
            double y = 0.4 + 0.2 * double(j) / SAMPLING;
            double z = 0.4 + 0.2 * double(k) / SAMPLING;
            coords.push_back(x);
            coords.push_back(y);
            coords.push_back(z);
          }
        }
      }
      break;
    }

    vertexIDs.resize(coords.size() / 3);
    participant.setMeshVertices(meshName, coords, vertexIDs);
    participant.initialize();
    dt = participant.getMaxTimeStepSize();

    // Run a step and read data expected to be f(x) = ax + by + cz + d
    BOOST_TEST(participant.isCouplingOngoing(), "Receiving participant must advance once.");

    participant.advance(dt);
    BOOST_TEST(!participant.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(vertexIDs.size());
    for (int i = 0; i < expected.size(); ++i) {
      expected(i) = a * coords[3 * i] + b * coords[3 * i + 1] + c * coords[3 * i + 2] + d;
    }
    Eigen::VectorXd readData(vertexIDs.size());

    dt = participant.getMaxTimeStepSize();
    participant.readData(meshName, dataName, vertexIDs, dt, readData);
    BOOST_CHECK(equals(expected, readData));
    participant.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
