#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(MultiCoupling)
/**
 * @brief Test to run a multi coupling with first order time interpolation.
 *
 * Provides a dt argument to the read function. A first order waveform is used.
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank))
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSamplingFirst)
{
  PRECICE_TEST();

  Participant precice(context.name, context.config(), 0, 1);

  std::string meshName, writeDataName;

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
  int    nSamples        = 4;
  int    iterations      = 0;
  double time            = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  double maxDt        = precice.getMaxTimeStepSize();
  double windowDt     = maxDt;
  double dt           = maxDt; // time step size desired by solver
  double currentDt    = dt;    // time step size used by solver
  double sampleDts[4] = {0.0, dt / 4.0, dt / 2.0, 3.0 * dt / 4.0};
  double readTime; // time where we are reading
  double sampleDt; // dt relative to timestep start, where we are sampling

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      windowStartTime = time;
      windowStartStep = timestep;
    }

    for (int j = 0; j < nSamples; j++) {
      sampleDt = sampleDts[j];
      readTime = time + sampleDt;

      for (auto &readDataPair : readDataPairs) {
        auto readDataName = readDataPair.first;
        auto readFunction = readDataPair.second;

        precice.readData(meshName, readDataName, {&vertexID, 1}, sampleDt, {&readData, 1});

        if (iterations == 0) { // always use constant extrapolation in first iteration (from initialize or writeData of second participant at end previous window).
          BOOST_TEST(readData == readFunction(time));
        } else if (iterations > 0) { // use linear interpolation in later iterations (additionally available writeData of second participant at end of this window).
          BOOST_TEST(readData == readFunction(readTime));
        } else {
          BOOST_TEST(false); // unreachable!
        }
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    maxDt     = precice.getMaxTimeStepSize();
    currentDt = dt > maxDt ? maxDt : dt;
    BOOST_CHECK(currentDt == windowDt); // no subcycling.
    time += currentDt;
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(maxDt);
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
  BOOST_TEST(timestep == nWindows);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // MultiCoupling

#endif // PRECICE_NO_MPI
