#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(ParallelCoupling)

/**
 * @brief Test to run a simple coupling where the number of time steps changes with the waveform iteration.
 *
 * Provides a dt argument to the read function. A first order waveform is used.
 */
BOOST_AUTO_TEST_CASE(WaveformSubcyclingWithDifferentNumberOfDts)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  typedef double (*DataFunction)(double);

  // The functions that are exchanged in iteration 1 and 2 respectively.
  DataFunction funcIterOne = [](double t) -> double {
    return (double) (1 + 2 * t) * (2 + 3 * t) * (3 + 4 * t);
  };
  DataFunction funcIterTwo = [](double t) -> double {
    return (double) (6 + 5 * t) * (1 + 6 * t) * (1 + 6 * t);
  };

  BOOST_TEST(funcIterOne(0) == funcIterTwo(0)); // because the data at WINDOW_START does not change over iterations.

  int         nTimeSteps;
  std::string meshName, writeDataName, readDataName;
  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "DataOne";
    readDataName  = "DataTwo";
    nTimeSteps    = 6;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "DataTwo";
    readDataName  = "DataOne";
    nTimeSteps    = 4;
  }

  double   writeData = 0;
  double   readData  = 0;
  double   v0[]      = {0, 0, 0};
  VertexID vertexID  = precice.setMeshVertex(meshName, v0);
  double   time      = 0;

  if (precice.requiresInitialData()) {
    writeData = funcIterOne(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  double maxDt     = precice.getMaxTimeStepSize();
  double windowDt  = maxDt;
  double dt        = windowDt / nTimeSteps;
  double currentDt = dt > maxDt ? maxDt : dt;

  double windowStartTime = 0;
  int    iterations      = 0;

  while (precice.isCouplingOngoing()) {

    if (precice.requiresWritingCheckpoint()) {
      windowStartTime = time;
      iterations      = 0;
    }

    if (context.isNamed("SolverOne")) {
      nTimeSteps = 6 + iterations; // increase the number of time steps by one for each iteration
    } else {
      BOOST_TEST(context.isNamed("SolverTwo"));
      nTimeSteps = 4 + iterations; // increase the number of time steps by one for each iteration
    }

    dt        = windowDt / nTimeSteps;
    maxDt     = precice.getMaxTimeStepSize();
    currentDt = dt > maxDt ? maxDt : dt;

    precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt, {&readData, 1});

    if (iterations == 0) { // in the first iteration of each window, we only have one sample of data. Therefore constant interpolation
      BOOST_TEST(readData == funcIterOne(windowStartTime));
    } else if (iterations == 1) { // Check that we receive the first function in the second iteration
      BOOST_TEST(readData == funcIterOne(time + currentDt));
    } else { // Check that we receive the second function in the third iteration
      BOOST_TEST(readData == funcIterTwo(time + currentDt));
    }

    // solve usually goes here. Dummy solve: Just sampling the two functions.
    time += currentDt;
    if (iterations == 0) { // Sends the first function in the first iteration
      writeData = funcIterOne(time);
    } else { //Otherwise send the second function
      writeData = funcIterTwo(time);
    }

    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(currentDt);

    if (precice.requiresReadingCheckpoint()) {
      time = windowStartTime;
      iterations++;
    }
  }

  precice.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // ParallelCoupling

#endif // PRECICE_NO_MPI
