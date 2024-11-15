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
 * @brief Test to run a simple coupling with zeroth order waveform subcycling.
 *
 * Provides a dt argument to the read function, uses zeroth order waveform for SolverOne and a first order waveform for SolverTwo. See ReadWriteScalarDataWithWaveformSubcyclingZero and ReadWriteScalarDataWithWaveformSubcyclingFirst for details on the non-mixed cases and expected behavior.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingMixed)
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

  double   writeData, readData;
  double   v0[]     = {0, 0, 0};
  VertexID vertexID = precice.setMeshVertex(meshName, v0);

  int    nSubsteps          = 4; // perform subcycling on solvers. 4 steps happen in each window.
  int    nWindows           = 5; // perform 5 windows.
  int    timestep           = 0;
  int    timestepCheckpoint = 0;
  double time               = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  double maxDt    = precice.getMaxTimeStepSize();
  double windowDt = maxDt;
  double dt       = windowDt / nSubsteps; // time step size desired by solver. E.g. 4 steps  with size 1/4
  dt += windowDt / nSubsteps / nSubsteps; // increase time step size such that we get a non-matching subcycling. E.g. 3 step with size 5/16 and 1 step with size 1/16.
  double currentDt = dt;                  // time step size used by solver
  double timeCheckpoint{0.0};
  int    iterations = 0;

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      timeCheckpoint     = time;
      timestepCheckpoint = timestep;
      iterations         = 0;
    }
    double readTime;
    if (context.isNamed("SolverOne")) {
      // zeroth order corresponds to end of window
      readTime = timeCheckpoint + windowDt;
    } else {
      BOOST_TEST(context.isNamed("SolverTwo"));
      // first order reads from interpolant inside window
      readTime = time + currentDt;
    }

    precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt, {&readData, 1});

    if (iterations == 0) { // in the first iteration of each window, use data from previous window.
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else {
      BOOST_TEST(readData == readFunction(readTime));
    }

    precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt / 2, {&readData, 1});

    if (iterations == 0) { // in the first iteration of each window, use data from previous window.
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else {                              // in the following iterations, use data at the end of window.
      if (context.isNamed("SolverOne")) { // in the following iterations, use data at the end of window.
        BOOST_TEST(readData == readFunction(readTime));
      } else { // in the following iterations we have two samples of data. Therefore linear interpolation
        BOOST_TEST(context.isNamed("SolverTwo"));
        BOOST_TEST(readData == readFunction(readTime - currentDt / 2));
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    timestep++;
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(currentDt);
    maxDt = precice.getMaxTimeStepSize();
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
BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // ParallelCoupling

#endif // PRECICE_NO_MPI
