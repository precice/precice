#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to run a simple coupling with subcycling.
 *
 *  Ensures that each time step provides its own data, but preCICE will only exchange data at the end of the window.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithSubcycling)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  SolverInterface precice(context.name, context.config(), 0, 1); // serial coupling, SolverOne first

  MeshID meshID;
  DataID writeDataID;
  DataID readDataID;

  typedef double (*DataFunction)(double, int);

  DataFunction dataOneFunction = [](double t, int idx) -> double {
    return (double) (2 + t + idx);
  };
  DataFunction dataTwoFunction = [](double t, int idx) -> double {
    return (double) (10 + t + idx);
  };
  DataFunction writeFunction;
  DataFunction readFunction;

  if (context.isNamed("SolverOne")) {
    meshID        = precice.getMeshID("MeshOne");
    writeDataID   = precice.getDataID("DataOne", meshID);
    writeFunction = dataOneFunction;
    readDataID    = precice.getDataID("DataTwo", meshID);
    readFunction  = dataTwoFunction;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshID        = precice.getMeshID("MeshTwo");
    writeDataID   = precice.getDataID("DataTwo", meshID);
    writeFunction = dataTwoFunction;
    readDataID    = precice.getDataID("DataOne", meshID);
    readFunction  = dataOneFunction;
  }

  int n_vertices = 1;

  std::vector<VertexID> vertexIDs(n_vertices, 0);
  std::vector<double>   writeData(n_vertices, 0);
  std::vector<double>   readData(n_vertices, 0);
  double                oldWriteData, oldReadData;

  vertexIDs[0] = precice.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());

  int    nSubsteps       = 4; // perform subcycling on solvers. 4 steps happen in each window.
  int    nWindows        = 5; // perform 5 windows.
  double maxDt           = precice.initialize();
  double windowDt        = maxDt;
  int    timestep        = 0;
  int    timewindow      = 0;
  double startTime       = 0;
  double windowStartTime = 0;
  int    windowStartStep = 0;
  int    iterations      = 0;
  double dt              = windowDt / (nSubsteps - 0.5); // Timestep length desired by solver. E.g. 4 steps with size 4/7. Fourth step will be restricted to 2/7 via preCICE steering to fit into the window.
  double expectedDts[]   = {4.0 / 7.0, 4.0 / 7.0, 4.0 / 7.0, 2.0 / 7.0};
  double currentDt       = dt; // Timestep length used by solver
  double time            = timestep * dt;

  if (precice.isActionRequired(precice::constants::actionWriteInitialData())) {
    for (int i = 0; i < n_vertices; i++) {
      writeData[i] = writeFunction(time, i);
      precice.writeScalarData(writeDataID, vertexIDs[i], writeData[i]);
    }
    precice.markActionFulfilled(precice::constants::actionWriteInitialData());
  }

  precice.initializeData();

  while (precice.isCouplingOngoing()) {
    if (precice.isActionRequired(precice::constants::actionWriteIterationCheckpoint())) {
      windowStartTime = time;
      windowStartStep = timestep;
      precice.markActionFulfilled(precice::constants::actionWriteIterationCheckpoint());
    }

    if (timestep % nSubsteps == 0) {
      BOOST_TEST(precice.isReadDataAvailable());
    } else {
      BOOST_TEST(!precice.isReadDataAvailable());
    }

    BOOST_TEST(readData.size() == n_vertices);
    // @todo split in SolverOne and SolverTwo?
    for (int i = 0; i < n_vertices; i++) {
      oldReadData = readData[i];
      if (precice.isReadDataAvailable()) {
        precice.readScalarData(readDataID, vertexIDs[i], readData[i]);
      }
      if (context.isNamed("SolverOne") && iterations == 0 && timestep == 0) {                      // special situation: SolverOne in its very first time window, first iteration, first time step
        BOOST_TEST(readData[i] != oldReadData);                                                    // update from uninitialized to initial data.
        BOOST_TEST(readData[i] == readFunction(startTime, i));                                     // use initial data only.
      } else if (context.isNamed("SolverOne") && iterations == 0) {                                // special situation: SolverOne gets the old data its first iteration for all time windows.
        BOOST_TEST(readData[i] == oldReadData);                                                    // ensure that read data stays the same from one step to the next, if not a new window is entered
        BOOST_TEST(readData[i] == readFunction(startTime + (timewindow) *windowDt, i));            // data at end of window was written by other solver.
      } else if (context.isNamed("SolverOne") && iterations == 1 && timestep == windowStartStep) { // special situation: SolverOne in its second iteration, first timestep of window
        BOOST_TEST(readData[i] != oldReadData);                                                    // ensure that read data stays the same from one step to the next, if not a new window is entered
        BOOST_TEST(readData[i] == readFunction(startTime + (timewindow + 1) * windowDt, i));       // data at end of window was written by other solver.
      } else if (context.isNamed("SolverTwo") && iterations == 0 && timestep == 0) {               // special situation: SolverTwo in its very first time window, first iteration, first time step
        BOOST_TEST(readData[i] != oldReadData);                                                    // update from uninitialized to initial data.
        BOOST_TEST(readData[i] == readFunction(startTime + (timewindow + 1) * windowDt, i));       // data at end of window was written by other solver.
      } else if (precice.isTimeWindowComplete()) {                                                 // moving to next window
        BOOST_TEST(readData[i] != oldReadData);                                                    // ensure that read data changes from one step to the next, if a new window is entered
        BOOST_TEST(readData[i] == readFunction(startTime + (timewindow + 1) * windowDt, i));       // data at end of window was written by other solver.
      } else if (not precice.isTimeWindowComplete()) {                                             // still iterating in the same window
        BOOST_TEST(readData[i] == oldReadData);                                                    // ensure that read data stays the same from one step to the next, if not a new window is entered
        BOOST_TEST(readData[i] == readFunction(startTime + (timewindow + 1) * windowDt, i));       // data at end of window was written by other solver.
      } else {                                                                                     // we should not enter this branch, because this would skip all tests.
        BOOST_TEST(false);
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    BOOST_TEST(currentDt == expectedDts[timestep % nSubsteps]);
    time += currentDt;
    if (precice.isWriteDataRequired(currentDt)) {
      BOOST_TEST(writeData.size() == n_vertices);
      for (int i = 0; i < n_vertices; i++) {
        oldWriteData = writeData[i];
        writeData[i] = writeFunction(time, i);
        precice.writeScalarData(writeDataID, vertexIDs[i], writeData[i]);
      }
    }
    maxDt     = precice.advance(currentDt);
    currentDt = dt > maxDt ? maxDt : dt;
    timestep++;
    if (precice.isActionRequired(precice::constants::actionReadIterationCheckpoint())) { // at end of window and we have to repeat it.
      iterations++;
      timestep = windowStartStep;
      time     = windowStartTime;
      precice.markActionFulfilled(precice::constants::actionReadIterationCheckpoint()); // this test does not care about checkpointing, but we have to make the action
    }
    if (precice.isTimeWindowComplete()) {
      timewindow++;
      iterations = 0;
    }
  }

  precice.finalize();
  BOOST_TEST(timestep == nWindows * nSubsteps);
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
