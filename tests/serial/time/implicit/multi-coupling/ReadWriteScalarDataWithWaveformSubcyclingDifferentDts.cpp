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
 * @brief Test to run a multi coupling with zeroth order waveform subcycling. Uses different time step sizes for each solver.
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank))
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingDifferentDts)
{
  PRECICE_TEST();

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
    readDataPairs.emplace_back(dataTwoName, dataTwoFunction);
    auto dataThreeName = "DataThree";
    readDataPairs.emplace_back(dataThreeName, dataThreeFunction);
    nSubsteps = 4;
  } else if (context.isNamed("SolverTwo")) {
    meshName         = "MeshTwo";
    writeDataName    = "DataTwo";
    writeFunction    = dataTwoFunction;
    auto dataOneName = "DataOne";
    readDataPairs.emplace_back(dataOneName, dataOneFunction);
    nSubsteps = 3;
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    meshName         = "MeshThree";
    writeDataName    = "DataThree";
    writeFunction    = dataThreeFunction;
    auto dataOneName = "DataOne";
    readDataPairs.emplace_back(dataOneName, dataOneFunction);
    nSubsteps = 2;
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
        BOOST_TEST(readData == readFunction(time + currentDt));
      }

      precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt / 2, {&readData, 1});

      if (iterations == 0) { // in the first iteration of each window, use data from previous window.
        BOOST_TEST(readData == readFunction(windowStartTime));
      } else { // in the following iterations, use data at the end of window.
        BOOST_TEST(readData == readFunction(time + currentDt / 2));
      }

      precice.readData(meshName, readDataName, {&vertexID, 1}, 0, {&readData, 1});

      if (iterations == 0) { // in the first iteration of each window, use data from previous window.
        BOOST_TEST(readData == readFunction(windowStartTime));
      } else { // in the following iterations, use data at the end of window.
        BOOST_TEST(readData == readFunction(time));
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    maxDt     = precice.getMaxTimeStepSize();
    currentDt = dt > maxDt ? maxDt : dt;
    BOOST_CHECK(math::equals(currentDt, windowDt / nSubsteps)); // subcycling.
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
