#ifndef PRECICE_NO_MPI

#include <precice/SolverInterface.hpp>
#include <vector>
#include "testing/Testing.hpp"
#include "utils/IntraComm.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(ParallelCubeConservative1To3)
{
  using precice::VertexID;
  using precice::testing::equals;

  PRECICE_TEST("SolverOneCubeConservative1To3"_on(1_rank), "SolverTwoCubeConservative1To3"_on(3_ranks).setupIntraComm());
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;
  double                dt;

  // Apply some forces (geometry described below)
  // They get spread to various ranks and tetra/triangle/edge
  double forceOnMidABC = 1.0; // Goes to rank2
  double forceOnMidACD = 0.5; // Goes to rank2
  //double forceOnA = 10; // Could go to any rank's duplicate of A, check the sum is conserved

  if (context.isNamed("SolverOneCubeConservative1To3")) {
    auto meshID = interface.getMeshID("MeshOne");
    auto dataID = interface.getDataID("DataOne", meshID);

    std::vector<double> coords;

    coords = {2. / 3, 1. / 3, 0,
              1. / 3, 2. / 3, 0.};

    vertexIDs.resize(coords.size() / 3);
    interface.setMeshVertices(meshID, vertexIDs.size(), coords.data(), vertexIDs.data());

    dt = interface.initialize();

    // Run a step and write forces
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    std::vector<double> values;
    values = {forceOnMidABC,
              forceOnMidACD};

    interface.writeBlockScalarData(dataID, values.size(), vertexIDs.data(), values.data());

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();
  } else { // SolverTwoCubeConservative1To3
    auto meshID = interface.getMeshID("MeshTwo");
    auto dataID = interface.getDataID("DataOne", meshID);

    std::vector<double> coords;

    /* 
    Let the cube be ABCDEFGH with A in 0,0,0, B in 1,0,0, C in 1,1,0, D in 0,1,0 
    and the rest similar but with z=1.
    Each rank consists of a pentahedron made of 2 tetra each. All contain A and G
    Rank 0 owns ADEGH
    Rank 1 owns ABEFG
    Rank 2 owns ABCDH
    */

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

    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");
    /*
    // Check expected VS read
    Eigen::VectorXd expected(3);
    Eigen::VectorXd readData(3);
    if (context.rank == 0) {
      // Force applied on (0.3, 0.5) to spread with respect to barycentric coordinates
      expected << 0.5, 0.3, 0.2;
    } else {
      // FOrce applied on (0.9, 0.2)
      expected << 0.1, 0.2, 0.7;
    }*/

    Eigen::VectorXd readData(5);
    interface.readBlockScalarData(dataID, vertexIDs.size(), vertexIDs.data(), readData.data());

    // map to global coordinates
    std::array<double, 8> forces{0, 0, 0, 0, 0, 0, 0, 0};
    std::array<double, 8> totalForces{0, 0, 0, 0, 0, 0, 0, 0};

    switch (context.rank) {
    case 0:
      forces[0] = readData[vertexIDs[0]]; //A
      forces[4] = readData[vertexIDs[1]]; //E
      forces[7] = readData[vertexIDs[2]]; //H
      forces[6] = readData[vertexIDs[3]]; //G
      forces[3] = readData[vertexIDs[4]]; //D
      break;
    case 1:
      forces[0] = readData[vertexIDs[0]]; //A
      forces[4] = readData[vertexIDs[1]]; //E
      forces[5] = readData[vertexIDs[2]]; //F
      forces[6] = readData[vertexIDs[3]]; //G
      forces[1] = readData[vertexIDs[4]]; //B
      break;
    case 2:
      forces[0] = readData[vertexIDs[0]]; //A
      forces[3] = readData[vertexIDs[1]]; //D
      forces[6] = readData[vertexIDs[2]]; //G
      forces[2] = readData[vertexIDs[3]]; //C
      forces[1] = readData[vertexIDs[4]]; //B
      break;

    default:
      break;
    }

    precice::utils::IntraComm::allreduceSum(forces, totalForces);

    BOOST_CHECK(equals(totalForces[0], 0.5));
    BOOST_CHECK(equals(totalForces[1], 0.5));
    BOOST_CHECK(equals(totalForces[2], 0.5));
    BOOST_CHECK(equals(totalForces[3], 0.5));
    BOOST_CHECK(equals(totalForces[4], 0.5));
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
