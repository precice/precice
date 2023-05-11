#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>
#include "precice/impl/ParticipantImpl.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(ParallelCubeConservative3To1)
{
  using precice::VertexID;
  using precice::testing::equals;

  PRECICE_TEST("SolverOneCubeConservative3To1"_on(3_ranks), "SolverTwoCubeConservative3To1"_on(1_rank));
  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  // Apply some forces (geometry described below)
  // They get spread to various ranks and tetra/triangle/edge
  // Each rank sends 2 of these
  double forceOnMidABC         = 1.0;
  double forceOnMidACD         = 0.5;
  double unbalancedForceOnGH   = 2.0; // 25% on G, 75% on H
  double forceOnMidAEGH        = 3.0;
  double forceNearC            = 7.0;
  double unbalancedForceOnAEGH = 7.0; // Distribution: 10%, 20%, 30%, 40

  std::vector<VertexID> vertexIDs;
  double                dt;

  if (context.isNamed("SolverOneCubeConservative3To1")) {
    auto meshName = "MeshOne";
    auto dataName = "DataOne";

    std::vector<double> coords;
    std::vector<double> values;

    // Each rank sends some "forces" on one point
    switch (context.rank) {
    case 0:
      coords = {2. / 3, 1. / 3, 0,
                1. / 3, 2. / 3, 0};
      values = {forceOnMidABC, forceOnMidACD};
      break;
    case 1:
      coords = {0.75, 1, 1,
                0.25, 0.5, 0.75};
      values = {unbalancedForceOnGH,
                forceOnMidAEGH};
      break;
    case 2:
      coords = {1.01, 1.01, 0.0,
                0.3, 0.7, 0.9};
      values = {forceNearC,
                unbalancedForceOnAEGH};
      break;
    default:
      break;
    }
    vertexIDs.resize(coords.size() / 3);
    interface.setMeshVertices(meshName, coords, vertexIDs);

    interface.initialize();
    dt = interface.getMaxTimeStepSize();

    // Run a step and write forces
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant must advance once.");

    interface.writeData(meshName, dataName, vertexIDs, values);

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant must advance only once.");
    interface.finalize();
  } else { // SolverTwoCubeConservative3To1
    auto meshName = "MeshTwo";
    auto dataName = "DataOne";

    std::vector<double> coords;

    // Unit cube made of 6 tetra, with 8 points
    coords = {0, 0, 0,
              1, 0, 0,
              1, 1, 0,
              0, 1, 0,
              0, 0, 1,
              1, 0, 1,
              1, 1, 1,
              0, 1, 1};

    vertexIDs.resize(coords.size() / 3);
    interface.setMeshVertices(meshName, coords, vertexIDs);

    VertexID v000 = vertexIDs[0];

    VertexID v100 = vertexIDs[1];
    VertexID v110 = vertexIDs[2];
    VertexID v010 = vertexIDs[3];
    VertexID v001 = vertexIDs[4];
    VertexID v101 = vertexIDs[5];
    VertexID v111 = vertexIDs[6];
    VertexID v011 = vertexIDs[7];

    interface.setMeshTetrahedron(meshName, v000, v001, v011, v111);
    interface.setMeshTetrahedron(meshName, v000, v010, v011, v111);
    interface.setMeshTetrahedron(meshName, v000, v001, v101, v111);
    interface.setMeshTetrahedron(meshName, v000, v100, v101, v111);
    interface.setMeshTetrahedron(meshName, v000, v010, v110, v111);
    interface.setMeshTetrahedron(meshName, v000, v100, v110, v111);

    auto &mesh = precice::testing::WhiteboxAccessor::impl(interface).mesh("MeshTwo");
    BOOST_REQUIRE(mesh.vertices().size() == 8);
    BOOST_REQUIRE(mesh.tetrahedra().size() == 6);
    interface.initialize();
    dt = interface.getMaxTimeStepSize();

    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant must advance once.");

    interface.advance(dt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant must advance only once.");

    Eigen::VectorXd readData(8);

    dt = interface.getMaxTimeStepSize();
    interface.readData(meshName, dataName, vertexIDs, dt, readData);
    BOOST_CHECK(equals(readData[0], forceOnMidABC / 3 + forceOnMidACD / 3 + forceOnMidAEGH / 4 + 0.1 * unbalancedForceOnAEGH));
    BOOST_CHECK(equals(readData[1], forceOnMidABC / 3));
    BOOST_CHECK(equals(readData[2], forceOnMidABC / 3 + forceOnMidACD / 3 + forceNearC));
    BOOST_CHECK(equals(readData[3], forceOnMidACD / 3));
    BOOST_CHECK(equals(readData[4], forceOnMidAEGH / 4 + 0.2 * unbalancedForceOnAEGH));
    BOOST_CHECK(equals(readData[5], 0.0));
    BOOST_CHECK(equals(readData[6], 0.75 * unbalancedForceOnGH + forceOnMidAEGH / 4 + 0.3 * unbalancedForceOnAEGH));
    BOOST_CHECK(equals(readData[7], 0.25 * unbalancedForceOnGH + forceOnMidAEGH / 4 + 0.4 * unbalancedForceOnAEGH));
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
