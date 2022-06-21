#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include "precice/impl/SolverInterfaceImpl.hpp"
#include <vector>
#include <iostream>
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(ParallelCube1To3)
{
  using precice::VertexID;
  using precice::testing::equals;

  PRECICE_TEST("SolverOne1_3"_on(1_rank), "SolverTwo1_3"_on(3_ranks));

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;
  double                a = 1, b = 2, c = 5, d = 1;

  if (context.isNamed("SolverOne1_3")) {
    auto meshID = interface.getMeshID("MeshOne");
    auto dataID = interface.getDataID("DataOne", meshID);

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
  } else { // SolverTwo1_3
    auto meshID = interface.getMeshID("MeshTwo");
    auto dataID = interface.getDataID("DataOne", meshID);

    std::vector<double> coords;

    // TODO: add more points

    switch (context.rank)
    {
    case 0:
      coords = {0.3, 0.3, 0.3,
                0.1, 0.1, 0.1};
      break;
    case 1:
      coords = {0.7, 0.7, 0.7,
                0.8, 0.9, 0.6};
      break;
    case 2:
      coords = {0.7, 0.3, 0.3};
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
      expected(i) = a * coords[3*i] + b * coords[3*i+1] + c * coords[3*i+2] + d;
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
