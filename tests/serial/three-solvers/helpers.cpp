#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "precice/SolverInterface.hpp"

/**
 * @brief Three solvers are coupled in a fork S2 <-> S1 <-> S3.
 *
 * Both couplings are explicit, solver 1 provides the mesh to the other two
 * solvers.
 */
void runTestThreeSolvers(std::string const &config, std::vector<int> expectedCallsOfAdvance, TestContext const &context)
{

  int callsOfAdvance = 0;

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface precice(context.name, config, 0, 1);

    auto meshAID = "MeshA";
    auto meshBID = "MeshB";
    precice.setMeshVertex(meshAID, Eigen::Vector2d(0, 0).data());
    precice.setMeshVertex(meshBID, Eigen::Vector2d(1, 1).data());

    if (precice.requiresInitialData()) {
    }
    double dt = precice.initialize();

    while (precice.isCouplingOngoing()) {
      if (precice.requiresWritingCheckpoint()) {
      }
      dt = precice.advance(dt);
      if (precice.requiresReadingCheckpoint()) {
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(0));
  } else if (context.isNamed("SolverTwo")) {
    SolverInterface precice(context.name, config, 0, 1);

    auto meshID = "MeshC";
    precice.setMeshVertex(meshID, Eigen::Vector2d(0, 0).data());

    if (precice.requiresInitialData()) {
    }
    double dt = precice.initialize();

    while (precice.isCouplingOngoing()) {
      if (precice.requiresWritingCheckpoint()) {
      }
      dt = precice.advance(dt);
      if (precice.requiresReadingCheckpoint()) {
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(1));
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    SolverInterface precice(context.name, config, 0, 1);

    auto meshID = "MeshD";
    precice.setMeshVertex(meshID, Eigen::Vector2d(0, 0).data());

    if (precice.requiresInitialData()) {
    }
    double dt = precice.initialize();

    while (precice.isCouplingOngoing()) {
      if (precice.requiresWritingCheckpoint()) {
      }
      dt = precice.advance(dt);
      if (precice.requiresReadingCheckpoint()) {
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(2));
  }
}

#endif
