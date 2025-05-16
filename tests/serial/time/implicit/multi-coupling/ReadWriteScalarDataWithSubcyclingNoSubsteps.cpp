#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(MultiCoupling)
/**
 * @brief Test to run a multi coupling with subcycling.
 *
 * Deactivates exchange of substeps.
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank))
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithSubcyclingNoSubsteps)
{
  PRECICE_TEST();

  Participant precice(context.name, context.config(), 0, 1);

  typedef double (*DataFunction)(double);
  std::vector<std::pair<std::string, DataFunction>> readDataPairs;

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 + t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (10 + t);
  };
  DataFunction dataThreeFunction = [](double t) -> double {
    return (double) (300 + t);
  };
  DataFunction writeFunction;

  int nSubsteps = 4;

  std::string meshName, writeDataName;
  if (context.isNamed("SolverOne")) {
    meshName         = "MeshOne";
    writeDataName    = "DataOne";
    writeFunction    = dataOneFunction;
    auto dataTwoName = "DataTwo";
    readDataPairs.emplace_back(dataTwoName, dataTwoFunction);
    auto dataThreeName = "DataThree";
    readDataPairs.emplace_back(dataThreeName, dataThreeFunction);
  } else if (context.isNamed("SolverTwo")) {
    meshName         = "MeshTwo";
    writeDataName    = "DataTwo";
    writeFunction    = dataTwoFunction;
    auto dataOneName = "DataOne";
    readDataPairs.emplace_back(dataOneName, dataOneFunction);
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    meshName         = "MeshThree";
    writeDataName    = "DataThree";
    writeFunction    = dataThreeFunction;
    auto dataOneName = "DataOne";
    readDataPairs.emplace_back(dataOneName, dataOneFunction);
  }

  double   writeData = 0;
  double   readData  = 0;
  double   v0[]      = {0, 0, 0};
  VertexID vertexID  = precice.setMeshVertex(meshName, v0);

  int    nWindows        = 5; // perform 5 windows.
  int    timestep        = 0;
  int    timewindow      = 0;
  double windowStartTime = 0;
  int    windowStartStep = 0;
  int    iterations      = 0;
  double time            = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  double windowDt = precice.getMaxTimeStepSize();
  double solverDt = windowDt / nSubsteps; // time step size desired by solver. E.g. 2 steps with size 1/2

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      windowStartTime = time;
      windowStartStep = timestep;
    }

    double preciceDt = precice.getMaxTimeStepSize();
    double currentDt = solverDt > preciceDt ? preciceDt : solverDt; // determine actual time step size; must fit into remaining time in window
    BOOST_CHECK(currentDt == windowDt / nSubsteps);                 // no subcycling.

    for (auto &readDataPair : readDataPairs) {
      auto readDataName = readDataPair.first;
      auto readFunction = readDataPair.second;

      precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt, {&readData, 1});
      if (iterations == 0 && timestep == 0) {                        // special situation: Both solvers are in their very first time windows, first iteration, first time step
        BOOST_TEST(readData == readFunction(0));                     // use initial data only.
      } else if (iterations == 0) {                                  // special situation: Both solvers get the old data for all time windows.
        BOOST_TEST(readData == readFunction(timewindow * windowDt)); // data at end of window was written by other solver.
      } else if (iterations > 0) {
        BOOST_TEST(readData == readFunction((timewindow + 1) * windowDt));
      } else { // we should not enter this branch, because this would skip all tests.
        BOOST_TEST(false);
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
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
      timewindow++;
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
BOOST_AUTO_TEST_SUITE_END() // MultiCoupling

#endif // PRECICE_NO_MPI
