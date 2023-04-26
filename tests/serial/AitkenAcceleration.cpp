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

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  Vector2d                 vertex{0.0, 0.0};

  if (context.isNamed("A")) {
    auto meshName = "A-Mesh";
    int  vertexID = interface.setMeshVertex(meshName, vertex.data());
    auto dataName = "Data";

    interface.initialize();
    double dt    = interface.getMaxTimeStepSize();
    double value = 1.0;
    interface.writeScalarData(meshName, dataName, vertexID, value);

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
    auto meshName = "B-Mesh";
    int  vertexID = interface.setMeshVertex(meshName, vertex.data());
    auto dataName = "Data";

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    interface.requiresWritingCheckpoint();
    interface.advance(dt);
    interface.requiresReadingCheckpoint();

    double value = -1.0;

    dt = interface.getMaxTimeStepSize();
    interface.readScalarData(meshName, dataName, vertexID, dt, value);
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
