#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(AitkenAcceleration)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  Vector2d                 vertex{0.0, 0.0};

  if (context.isNamed("A")) {
    const precice::MeshID meshID   = interface.getMeshID("A-Mesh");
    int                   vertexID = interface.setMeshVertex(meshID, vertex.data());
    int                   dataID   = interface.getDataID("Data", meshID);

    double dt    = interface.initialize();
    double value = 1.0;
    interface.writeScalarData(dataID, vertexID, value);

    interface.requiresWritingCheckpoint();
    interface.advance(dt);
    interface.requiresReadingCheckpoint();

    interface.requiresWritingCheckpoint();
    interface.advance(dt);
    interface.requiresReadingCheckpoint();

    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    const precice::MeshID meshID   = interface.getMeshID("B-Mesh");
    int                   vertexID = interface.setMeshVertex(meshID, vertex.data());
    int                   dataID   = interface.getDataID("Data", meshID);

    double dt = interface.initialize();
    interface.requiresWritingCheckpoint();
    interface.advance(dt);
    interface.requiresReadingCheckpoint();

    double value = -1.0;
    interface.readScalarData(dataID, vertexID, value);
    BOOST_TEST(value == 0.1); // due to initial underrelaxation

    interface.requiresWritingCheckpoint();
    interface.advance(dt);
    interface.requiresReadingCheckpoint();

    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
