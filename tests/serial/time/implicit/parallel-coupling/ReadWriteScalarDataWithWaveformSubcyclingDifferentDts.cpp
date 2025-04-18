#ifndef PRECICE_NO_MPI

#include "math/differences.hpp"
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
 * @brief Test to run a simple coupling with zeroth order waveform subcycling. Uses different time step sizes for both solvers.
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingDifferentDts)
{
  PRECICE_TEST();

  Participant precice(context.name, context.config(), 0, 1);

  std::string meshName;
  std::string writeDataName;
  std::string readDataName;

  typedef double (*DataFunction)(double);

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 + t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (10 + t);
  };
  DataFunction writeFunction;
  DataFunction readFunction;

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

  int nSubsteps; // perform subcycling on solvers. nSubsteps steps happen in each window.
  if (context.isNamed("SolverOne")) {
    nSubsteps = 4;
  } else {
    nSubsteps = 3;
  }
  int    nWindows = 5; // perform 5 windows.
  int    timestep = 0;
  double time     = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  double maxDt              = precice.getMaxTimeStepSize();
  double windowDt           = maxDt;
  int    timestepCheckpoint = 0;
  double dt                 = windowDt / nSubsteps; // Timestep length desired by solver. E.g. 4 steps  with size 1/4, 3 steps with size 1/3
  dt += windowDt / nSubsteps / nSubsteps;           // increase timestep such that we get a non-matching subcycling. E.g. 3 steps with size 5/16 and 1 steps with size 1/16; 2 steps with size 4/9 and 1 step with size 1/9
  double currentDt      = dt;                       // Timestep length used by solver
  double timeCheckpoint = 0.0;
  int    iterations     = 0;

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      timeCheckpoint     = time;
      timestepCheckpoint = timestep;
      iterations         = 0;
    }
    precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt, {&readData, 1});

    if (iterations == 0) { // in the first iteration of each window, use data from previous window.
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else { // in the following iterations, use data at the end of window.
      BOOST_TEST(readData == readFunction(time + currentDt));
    }

    precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt / 2, {&readData, 1});

    if (iterations == 0) { // in the first iteration of each window, use data from previous window.
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else { // in the following iterations, use data at the end of window.
      BOOST_TEST(readData == readFunction(time + currentDt / 2));
    }

    precice.readData(meshName, readDataName, {&vertexID, 1}, 0, {&readData, 1});

    if (iterations == 0) { // in the first iteration of each window, use data from previous window.
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else { // in the following iterations, use data at the end of window.
      BOOST_TEST(readData == readFunction(time));
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    timestep++;
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(currentDt);
    double maxDt = precice.getMaxTimeStepSize();
    if (precice.requiresReadingCheckpoint()) {
      time     = timeCheckpoint;
      timestep = timestepCheckpoint;
      iterations++;
    }
    currentDt = dt > maxDt ? maxDt : dt;
  }

  precice.finalize();
  BOOST_TEST(timestep == nWindows * nSubsteps);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
