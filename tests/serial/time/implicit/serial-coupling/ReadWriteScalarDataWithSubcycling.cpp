#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to run a simple coupling with subcycling.
 *
 * Ensures that each time step provides its own data, but preCICE only exchanges data at the end of the window.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithSubcycling)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  SolverInterface precice(context.name, context.config(), 0, 1);

  typedef double (*DataFunction)(double);

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 + t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (10 + t);
  };
  DataFunction writeFunction;
  DataFunction readFunction;

  std::string meshName, writeDataID, readDataID;
  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataID   = "DataOne"; //  meshName
    writeFunction = dataOneFunction;
    readDataID    = "DataTwo"; //  meshName
    readFunction  = dataTwoFunction;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataID   = "DataTwo"; //  meshName
    writeFunction = dataTwoFunction;
    readDataID    = "DataOne"; //  meshName
    readFunction  = dataOneFunction;
  }

  double   writeData = 0;
  double   readData  = 0;
  VertexID vertexID  = precice.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());

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
    precice.writeScalarData(meshName, writeDataID, vertexID, writeData);
  }

  double maxDt = precice.initialize();
  BOOST_TEST(maxDt == 2.0); // use window size != 1.0 to be able to detect more possible bugs
  double windowDt      = maxDt;
  double dt            = windowDt / (nSubsteps - 0.5);                 // Solver always tries to do a timestep of fixed size.
  double expectedDts[] = {4.0 / 7.0, 4.0 / 7.0, 4.0 / 7.0, 2.0 / 7.0}; // If solver uses timestep size of 4/7, fourth step will be restricted to 2/7 via preCICE steering to fit into the window.
  double currentDt     = dt > maxDt ? maxDt : dt;                      // determine actual timestep length; must fit into remaining time in window

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      windowStartTime = time;
      windowStartStep = timestep;
    }

    precice.readScalarData(meshName, readDataID, vertexID, readData);

    if (context.isNamed("SolverOne") && iterations == 0) {                     // special situation for serial coupling: SolverOne gets the old data in its first iteration for all time windows.
      BOOST_TEST(readData == readFunction(startTime + timewindow * windowDt)); // zeroth window: Initial Data from SolverTwo; following windows: data at end of window was written by SolverTwo.
    } else {
      BOOST_TEST(readData == readFunction(startTime + (timewindow + 1) * windowDt));
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    BOOST_TEST(currentDt == expectedDts[timestep % nSubsteps]);
    time += currentDt;
    writeData = writeFunction(time);
    precice.writeScalarData(meshName, writeDataID, vertexID, writeData);
    maxDt     = precice.advance(currentDt);
    currentDt = dt > maxDt ? maxDt : dt;
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
