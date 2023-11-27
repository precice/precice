#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include <fstream>
#include "precice/precice.hpp"

void subcyclingWithNSteps(TestContext const &context, int nSubsteps, std::vector<int> expectedSteps, bool useAdvancedDtStrategy)
{
  Participant precice(context.name, context.config(), 0, 1);

  std::string meshName, writeDataName, readDataName;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "DataOne";
    readDataName  = "DataTwo";
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "DataTwo";
    readDataName  = "DataOne";
  }

  double writeData, readData;

  double   v0[]     = {0, 0, 0};
  VertexID vertexID = precice.setMeshVertex(meshName, v0);

  int    timestep   = 0;
  int    timewindow = 0;
  double time       = 0;

  if (precice.requiresInitialData()) {
    writeData = 1; // don't care
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  BOOST_TEST(precice.getMaxTimeStepSize() == 0.2);
  double windowDt = precice.getMaxTimeStepSize();
  double solverDt = windowDt / nSubsteps;
  int    didSteps = 0;
  int    nWindows = 0;

  while (precice.isCouplingOngoing()) {
    if (precice.isTimeWindowComplete()) {
      BOOST_TEST(didSteps == expectedSteps[nWindows]);
      didSteps = 0; // reset counter for next window
      nWindows++;
    }
    double preciceDt = precice.getMaxTimeStepSize();

    // Correct strategy to compute solver dt that users should apply to avoid PRECICE_ERROR
    double currentDt;

    if (not useAdvancedDtStrategy) {
      // Simple strategy for determining time step size fails for many substeps
      currentDt = solverDt > preciceDt ? preciceDt : solverDt;
    } else {
      // Advanced strategy for determining time step size considers round of errors
      double tol = 100 * math::NUMERICAL_ZERO_DIFFERENCE;

      if (abs(preciceDt - solverDt) < tol) {
        currentDt = preciceDt;
      } else {
        currentDt = solverDt > preciceDt ? preciceDt : solverDt;
      }
      currentDt = solverDt > preciceDt ? preciceDt : solverDt;
    }

    precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt, {&readData, 1});
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(currentDt);
    didSteps++;
  }

  precice.finalize();
}

#endif
