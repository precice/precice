#ifndef PRECICE_NO_MPI

#include "math/differences.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(ParallelCoupling)

/**
 * @brief Test to run a simple coupling with subcycling.
 *
 * This is a smoke test to reproduce the scenario explained in https://github.com/precice/precice/issues/1866
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithSubcyclingSmallSteps)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  std::string meshName, writeDataName, readDataName;
  int         nSubsteps;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "DataOne";
    readDataName  = "DataTwo";
    nSubsteps     = 6400; // We get problems now for 6400
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "DataTwo";
    readDataName  = "DataOne";
    nSubsteps     = 1;
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
  double           windowDt = precice.getMaxTimeStepSize();
  double           solverDt = windowDt / nSubsteps;
  std::vector<int> expectedSteps{nSubsteps, nSubsteps, nSubsteps, nSubsteps, nSubsteps};
  int              didSteps = 0;
  int              nWindows = 0;

  while (precice.isCouplingOngoing()) {
    if (precice.isTimeWindowComplete()) {
      BOOST_TEST(didSteps == expectedSteps[nWindows]);
      didSteps = 0;
      nWindows++;
    }
    double preciceDt = precice.getMaxTimeStepSize();

    // Correct strategy to compute solver dt that users should apply to avoid PRECICE_ERROR
    double tol = math::NUMERICAL_ZERO_DIFFERENCE;
    double currentDt;
    if (abs(preciceDt - solverDt) < tol) {
      currentDt = preciceDt;
    } else {
      currentDt = solverDt > preciceDt ? preciceDt : solverDt;
    }
    currentDt = solverDt > preciceDt ? preciceDt : solverDt; // Bypasses safeguard above and triggers a PRECICE_ERROR after some time windows. As soon as https://github.com/precice/precice/issues/280 is merged, we should check for the error.

    precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt, {&readData, 1});
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(currentDt);
    didSteps++;
  }

  precice.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
