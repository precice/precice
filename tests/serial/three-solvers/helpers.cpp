#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "precice/Participant.hpp"

/**
 * @brief Three solvers are coupled in a fork S2 <-> S1 <-> S3.
 *
 * Both couplings are explicit, solver 1 provides the mesh to the other two
 * solvers.
 */
void runTestThreeSolvers(std::string const &config, std::vector<int> expectedCallsOfAdvance, TestContext const &context)
{

  int callsOfAdvance = 0;

  double v0[] = {0, 0};
  double v1[] = {1, 1};

  if (context.isNamed("SolverOne")) {
    precice::Participant precice(context.name, config, 0, 1);

    auto meshAID = "MeshA";
    precice.setMeshVertex(meshAID, v0);

    if (precice.requiresInitialData()) {
    }
    precice.initialize();
    double dt = precice.getMaxTimeStepSize();

    while (precice.isCouplingOngoing()) {
      if (precice.requiresWritingCheckpoint()) {
      }
      precice.advance(dt);
      dt = precice.getMaxTimeStepSize();
      if (precice.requiresReadingCheckpoint()) {
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(0));
  } else if (context.isNamed("SolverTwo")) {
    Participant precice(context.name, config, 0, 1);

    auto meshName = "MeshC";
    precice.setMeshVertex(meshName, v0);

    if (precice.requiresInitialData()) {
    }
    precice.initialize();
    double dt = precice.getMaxTimeStepSize();

    while (precice.isCouplingOngoing()) {
      if (precice.requiresWritingCheckpoint()) {
      }
      precice.advance(dt);
      dt = precice.getMaxTimeStepSize();
      if (precice.requiresReadingCheckpoint()) {
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(1));
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    Participant precice(context.name, config, 0, 1);

    auto meshName = "MeshD";
    precice.setMeshVertex(meshName, v0);

    if (precice.requiresInitialData()) {
    }
    precice.initialize();
    double dt = precice.getMaxTimeStepSize();

    while (precice.isCouplingOngoing()) {
      if (precice.requiresWritingCheckpoint()) {
      }
      precice.advance(dt);
      dt = precice.getMaxTimeStepSize();
      if (precice.requiresReadingCheckpoint()) {
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(2));
  }
}

#endif
