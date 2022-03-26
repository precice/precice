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
  std::string writeIterCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterCheckpoint(constants::actionReadIterationCheckpoint());
  std::string writeInitData(constants::actionWriteInitialData());

  int callsOfAdvance = 0;

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface precice(context.name, config, 0, 1);

    int meshAID = precice.getMeshID("MeshA");
    int meshBID = precice.getMeshID("MeshB");
    precice.setMeshVertex(meshAID, Eigen::Vector2d(0, 0).data());
    precice.setMeshVertex(meshBID, Eigen::Vector2d(1, 1).data());
    double dt = precice.initialize();

    if (precice.isActionRequired(writeInitData)) {
      precice.markActionFulfilled(writeInitData);
    }
    precice.initializeData();

    while (precice.isCouplingOngoing()) {
      if (precice.isActionRequired(writeIterCheckpoint)) {
        precice.markActionFulfilled(writeIterCheckpoint);
      }
      dt = precice.advance(dt);
      if (precice.isActionRequired(readIterCheckpoint)) {
        precice.markActionFulfilled(readIterCheckpoint);
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(0));
  } else if (context.isNamed("SolverTwo")) {
    SolverInterface precice(context.name, config, 0, 1);

    MeshID meshID = precice.getMeshID("MeshC");
    precice.setMeshVertex(meshID, Eigen::Vector2d(0, 0).data());
    double dt = precice.initialize();

    if (precice.isActionRequired(writeInitData)) {
      precice.markActionFulfilled(writeInitData);
    }
    precice.initializeData();

    while (precice.isCouplingOngoing()) {
      if (precice.isActionRequired(writeIterCheckpoint)) {
        precice.markActionFulfilled(writeIterCheckpoint);
      }
      dt = precice.advance(dt);
      if (precice.isActionRequired(readIterCheckpoint)) {
        precice.markActionFulfilled(readIterCheckpoint);
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(1));
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    SolverInterface precice(context.name, config, 0, 1);

    MeshID meshID = precice.getMeshID("MeshD");
    precice.setMeshVertex(meshID, Eigen::Vector2d(0, 0).data());
    double dt = precice.initialize();

    if (precice.isActionRequired(writeInitData)) {
      precice.markActionFulfilled(writeInitData);
    }
    precice.initializeData();

    while (precice.isCouplingOngoing()) {
      if (precice.isActionRequired(writeIterCheckpoint)) {
        precice.markActionFulfilled(writeIterCheckpoint);
      }
      dt = precice.advance(dt);
      if (precice.isActionRequired(readIterCheckpoint)) {
        precice.markActionFulfilled(readIterCheckpoint);
      }
      callsOfAdvance++;
    }
    precice.finalize();
    BOOST_TEST(callsOfAdvance == expectedCallsOfAdvance.at(2));
  }
}

#endif