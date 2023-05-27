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

int solverOneNSubsteps   = 4;
int solverTwoNSubsteps   = 3;
int solverThreeNSubsteps = 2;

double matchTimeFromOtherSolver(double thisTime, double windowStartTime, int otherSubsteps, double *otherTimeGrid)
{
  double returnedTime = windowStartTime;

  if (math::equals(windowStartTime, thisTime)) { // we are at the beginning of the window
    return returnedTime;
  }

  for (int i = 0; i < otherSubsteps; i++) { // step through all times on the grid of the other solver
    double relativeDt = otherTimeGrid[i];
    returnedTime      = windowStartTime + relativeDt;                  // point in time on grid of other solver
    if (math::greaterEquals(windowStartTime + relativeDt, thisTime)) { // thisTime lies between windowStartTime+otherTimeGrid[i-1] and windowStartTime+otherTimeGrid[i]
      return returnedTime;                                             // return windowStartTime+otherTimeGrid[i]
    }
  }
  BOOST_TEST(false); // unreachable
  return -1;
}

// helper to map a time to the corresponding time on the time grid of solver one. Helps to determine the expected value in the constant interpolation
double solverOneTime(double time, double windowStartTime)
{
  BOOST_TEST(solverOneNSubsteps == 4);
  double relativeDts[4] = {0.5, 1.0, 1.5, 2.0};
  return matchTimeFromOtherSolver(time, windowStartTime, solverOneNSubsteps, relativeDts);
}

// helper to map a time to the corresponding time on the time grid of solver two. Helps to determine the expected value in the constant interpolation
double solverTwoTime(double time, double windowStartTime)
{
  BOOST_TEST(solverTwoNSubsteps == 3);
  double relativeDts[3] = {2.0 / 3, 4.0 / 3, 2.0};
  return matchTimeFromOtherSolver(time, windowStartTime, solverTwoNSubsteps, relativeDts);
}

// helper to map a time to the corresponding time on the time grid of solver three. Helps to determine the expected value in the constant interpolation
double solverThreeTime(double time, double windowStartTime)
{
  BOOST_TEST(solverThreeNSubsteps == 2);
  double relativeDts[2] = {1.0, 2.0};
  return matchTimeFromOtherSolver(time, windowStartTime, solverThreeNSubsteps, relativeDts);
}

/**
 * @brief Test to run a multi coupling with zeroth order waveform subcycling. Uses different time step sizes for each solver.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingDifferentDts)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  std::string meshName;
  std::string writeDataName;

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

  int nSubsteps; // let three solvers use different time step sizes

  if (context.isNamed("SolverOne")) {
    meshName         = "MeshOne";
    writeDataName    = "DataOne";
    writeFunction    = dataOneFunction;
    auto dataTwoName = "DataTwo";
    readDataPairs.push_back(std::make_pair(dataTwoName, dataTwoFunction));
    auto dataThreeName = "DataThree";
    readDataPairs.push_back(std::make_pair(dataThreeName, dataThreeFunction));
    nSubsteps = solverOneNSubsteps;
  } else if (context.isNamed("SolverTwo")) {
    meshName         = "MeshTwo";
    writeDataName    = "DataTwo";
    writeFunction    = dataTwoFunction;
    auto dataOneName = "DataOne";
    readDataPairs.push_back(std::make_pair(dataOneName, dataOneFunction));
    nSubsteps = solverTwoNSubsteps;
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    meshName         = "MeshThree";
    writeDataName    = "DataThree";
    writeFunction    = dataThreeFunction;
    auto dataOneName = "DataOne";
    readDataPairs.push_back(std::make_pair(dataOneName, dataOneFunction));
    nSubsteps = solverThreeNSubsteps;
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
  int    nSamples        = 4;
  int    iterations      = 0;
  double time            = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  double maxDt     = precice.getMaxTimeStepSize();
  double windowDt  = maxDt;
  double dt        = windowDt / nSubsteps; // Timestep length desired by solver. E.g. 2 steps with size 1/2
  double currentDt = dt;                   // Timestep length used by solver

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      windowStartTime = time;
      windowStartStep = timestep;
      iterations      = 0;
    }

    for (auto readDataPair : readDataPairs) {
      auto readDataName = readDataPair.first;
      auto readFunction = readDataPair.second;

      precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt, {&readData, 1});

      if (iterations == 0) { // in the first iteration of each window, use data from previous window.
        BOOST_TEST(readData == readFunction(windowStartTime));
      } else { // in the following iterations, use data at the end of window.
        double readTime;
        if (readDataName == "DataOne") {
          readTime = solverOneTime(time + currentDt, windowStartTime);
        } else if (readDataName == "DataTwo") {
          readTime = solverTwoTime(time + currentDt, windowStartTime);
        } else if (readDataName == "DataThree") {
          readTime = solverThreeTime(time + currentDt, windowStartTime);
        } else {
          BOOST_TEST(false);
        }
        BOOST_TEST(readData == readFunction(readTime));
      }

      precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt / 2, {&readData, 1});

      if (iterations == 0) { // in the first iteration of each window, use data from previous window.
        BOOST_TEST(readData == readFunction(windowStartTime));
      } else { // in the following iterations, use data at the end of window.
        double readTime;
        if (readDataName == "DataOne") {
          readTime = solverOneTime(time + currentDt / 2, windowStartTime);
        } else if (readDataName == "DataTwo") {
          readTime = solverTwoTime(time + currentDt / 2, windowStartTime);
        } else if (readDataName == "DataThree") {
          readTime = solverThreeTime(time + currentDt / 2, windowStartTime);
        } else {
          BOOST_TEST(false);
        }
        BOOST_TEST(readData == readFunction(readTime));
      }

      precice.readData(meshName, readDataName, {&vertexID, 1}, 0, {&readData, 1});

      if (iterations == 0) { // in the first iteration of each window, use data from previous window.
        BOOST_TEST(readData == readFunction(windowStartTime));
      } else { // in the following iterations, use data at the end of window.
        double readTime;
        if (readDataName == "DataOne") {
          readTime = solverOneTime(time, windowStartTime);
        } else if (readDataName == "DataTwo") {
          readTime = solverTwoTime(time, windowStartTime);
        } else if (readDataName == "DataThree") {
          readTime = solverThreeTime(time, windowStartTime);
        } else {
          BOOST_TEST(false);
        }
        BOOST_TEST(readData == readFunction(readTime));
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(currentDt);
    double maxDt = precice.getMaxTimeStepSize();
    currentDt    = dt > maxDt ? maxDt : dt;
    BOOST_CHECK(currentDt == windowDt / nSubsteps); // no subcycling.
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
