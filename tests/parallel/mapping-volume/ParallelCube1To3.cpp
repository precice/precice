#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <iostream>
#include <precice/SolverInterface.hpp>
#include <vector>
#include "precice/impl/SolverInterfaceImpl.hpp"
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(ParallelCube1To3)
{
  using precice::VertexID;
  using precice::testing::equals;

  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(3_ranks));

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;
  double                a = 1, b = 2, c = 5, d = 1;

  if (context.isNamed("SolverOne")) {
    auto meshID = "MeshOne";
    auto dataID = "DataOne"; //  meshID

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
    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());

    VertexID v000 = vertexIDs[0];
    VertexID v001 = vertexIDs[1];
    VertexID v010 = vertexIDs[2];
    VertexID v011 = vertexIDs[3];
    VertexID v100 = vertexIDs[4];
    VertexID v101 = vertexIDs[5];
    VertexID v110 = vertexIDs[6];
    VertexID v111 = vertexIDs[7];

    interface.setMeshTetrahedron(meshID, v000, v001, v011, v111);
    interface.setMeshTetrahedron(meshID, v000, v010, v011, v111);
    interface.setMeshTetrahedron(meshID, v000, v001, v101, v111);
    interface.setMeshTetrahedron(meshID, v000, v100, v101, v111);
    interface.setMeshTetrahedron(meshID, v000, v010, v110, v111);
    interface.setMeshTetrahedron(meshID, v000, v100, v110, v111);

    auto &mesh = precice::testing::WhiteboxAccessor::impl(interface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 8);
    BOOST_REQUIRE(mesh.tetrahedra().size() == 6);

    dt = interface.initialize();

    // Run a step and write data with f(x) = ax + by + cz + d
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    values = {d,
              c + d,
              b + d,
              b + c + d,
              a + d,
              a + c + d,
              a + b + d,
              c + a + b + d};

    interface.writeBlockScalarData(dataID, values.size(), vertexIDs.data(), values.data());

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();
  } else { // SolverTwo
    auto meshID = "MeshTwo";
    auto dataID = "DataOne"; //  meshID

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
    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());
    dt = interface.initialize();

    // Run a step and read data expected to be f(x) = ax + by + cz + d
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    // Check expected VS read
    Eigen::VectorXd expected(vertexIDs.size());
    for (int i = 0; i < expected.size(); ++i) {
      expected(i) = a * coords[3 * i] + b * coords[3 * i + 1] + c * coords[3 * i + 2] + d;
    }
    Eigen::VectorXd readData(vertexIDs.size());

    interface.readBlockScalarData(dataID, expected.size(), vertexIDs.data(), readData.data());
    BOOST_CHECK(equals(expected, readData));
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
