#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(ParallelCoupling)

/**
 * @brief Test to run a simple coupling with subcycling.
 *
 * Ensures that each time step provides its own data, but preCICE only exchanges data at the end of the window.
 *
 * Deactivates exchange of substeps. Uses first degree interpolation. Uses data initialization for both participants.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithSubcyclingNoSubstepsUseInitFirst)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  typedef double (*DataFunction)(double);

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 + t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (10 + t);
  };
  DataFunction writeFunction;
  DataFunction readFunction;

  std::string meshName, writeDataName, readDataName;
  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "DataOne";
    writeFunction = dataOneFunction;
    readDataName  = "DataTwo";
    readFunction  = dataTwoFunction;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "DataTwo";
    writeFunction = dataTwoFunction;
    readDataName  = "DataOne";
    readFunction  = dataOneFunction;
  }

  double   writeData = 0;
  double   readData  = 0;
  double   v0[]      = {0, 0, 0};
  VertexID vertexID  = precice.setMeshVertex(meshName, v0);

  int    nSubsteps       = 4; // perform subcycling on solvers. 4 steps happen in each window.
  int    nWindows        = 5; // perform 5 windows.
  int    timestep        = 0;
  int    timewindow      = 0;
  double startTime       = 0;
  double windowStartTime = 0;
  int    windowStartStep = 0;
  int    iterations      = 0;
  double time            = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  BOOST_TEST(precice.getMaxTimeStepSize() == 2.0); // use window size != 1.0 to be able to detect more possible bugs
  double windowDt      = precice.getMaxTimeStepSize();
  double solverDt      = windowDt / (nSubsteps - 0.5);                 // Solver always tries to do a timestep of fixed size.
  double expectedDts[] = {4.0 / 7.0, 4.0 / 7.0, 4.0 / 7.0, 2.0 / 7.0}; // If solver uses timestep size of 4/7, fourth step will be restricted to 2/7 via preCICE steering to fit into the window.

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      windowStartTime = time;
      windowStartStep = timestep;
    }
    double preciceDt = precice.getMaxTimeStepSize();
    double currentDt = solverDt > preciceDt ? preciceDt : solverDt; // determine actual time step size; must fit into remaining time in window

    precice.readData(meshName, readDataName, {&vertexID, 1}, 0 * currentDt, {&readData, 1});

    if (iterations == 0) {                                                     // special situation: Both solvers are in their very first time windows, first iteration, first time step
      BOOST_TEST(readData == readFunction(startTime + timewindow * windowDt)); // zeroth window: Initial Data from SolverTwo; following windows: data at end of window was written by SolverTwo.
    } else {
      BOOST_TEST(readData == readFunction(time + 0 * currentDt));
    }

    precice.readData(meshName, readDataName, {&vertexID, 1}, 0.5 * currentDt, {&readData, 1});

    if (iterations == 0) {                                                     // special situation: Both solvers are in their very first time windows, first iteration, first time step
      BOOST_TEST(readData == readFunction(startTime + timewindow * windowDt)); // zeroth window: Initial Data from SolverTwo; following windows: data at end of window was written by SolverTwo.
    } else {
      BOOST_TEST(readData == readFunction(time + 0.5 * currentDt));
    }

    precice.readData(meshName, readDataName, {&vertexID, 1}, 1.0 * currentDt, {&readData, 1});

    if (iterations == 0) {                                                     // special situation: Both solvers are in their very first time windows, first iteration, first time step
      BOOST_TEST(readData == readFunction(startTime + timewindow * windowDt)); // zeroth window: Initial Data from SolverTwo; following windows: data at end of window was written by SolverTwo.
    } else {
      BOOST_TEST(readData == readFunction(time + 1.0 * currentDt));
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    BOOST_TEST(currentDt == expectedDts[timestep % nSubsteps]);
    time += currentDt;
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(currentDt);
    timestep++;
    if (precice.requiresReadingCheckpoint()) { // at end of window and we have to repeat it.
      iterations++;
      timestep = windowStartStep;
      time     = windowStartTime;
    }
    if (precice.isTimeWindowComplete()) {
      iterations++;
      timewindow++;
      BOOST_TEST(iterations == 3);
      iterations = 0;
    }
  }

  precice.finalize();
  BOOST_TEST(timestep == nWindows * nSubsteps);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
