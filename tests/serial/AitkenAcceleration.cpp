#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(AitkenAcceleration)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;
  using namespace precice::constants;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  Vector2d                 vertex{0.0, 0.0};

  if (context.isNamed("A")) {
    const precice::MeshID meshID   = interface.getMeshID("A-Mesh");
    int                   vertexID = interface.setMeshVertex(meshID, vertex.data());
    int                   dataID   = interface.getDataID("Data", meshID);

    double dt    = interface.initialize();
    double value = 1.0;
    interface.writeScalarData(dataID, vertexID, value);

    interface.markActionFulfilled(actionWriteIterationCheckpoint());
    interface.advance(dt);
    interface.markActionFulfilled(actionReadIterationCheckpoint());
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    const precice::MeshID meshID   = interface.getMeshID("B-Mesh");
    int                   vertexID = interface.setMeshVertex(meshID, vertex.data());
    int                   dataID   = interface.getDataID("Data", meshID);

    double dt = interface.initialize();
    interface.markActionFulfilled(actionWriteIterationCheckpoint());
    interface.advance(dt);

    double value = -1.0;
    interface.readScalarData(dataID, vertexID, value);
    BOOST_TEST(value == 0.1); // due to initial underrelaxation

    interface.markActionFulfilled(actionReadIterationCheckpoint());
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
