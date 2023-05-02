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
 * @brief Test to run a simple coupling with first order waveform subcycling.
 *
 * Does not call initializeData and therefore automatically uses 0 initial data.
 * Provides a dt argument to the read function. A first order waveform is used.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSamplingFirstNoInit)
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
  VertexID vertexID  = precice.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());

  int nWindows = 5; // perform 5 windows.
  precice.initialize();
  double maxDt           = precice.getMaxTimeStepSize();
  double windowDt        = maxDt;
  int    timestep        = 0;
  int    timewindow      = 0;
  double windowStartTime = 0;
  int    windowStartStep = 0;
  double dt              = maxDt; // time step size desired by solver
  double currentDt       = dt;    // time step size used by solver
  double time            = timestep * dt;
  double sampleDts[4]    = {0.0, dt / 4.0, dt / 2.0, 3.0 * dt / 4.0};
  int    nSamples        = 4;
  int    iterations      = 0;
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
      precice.readScalarData(meshName, readDataName, vertexID, sampleDt, readData);

      if (context.isNamed("SolverOne") && iterations == 0 && timewindow == 0) { // use zero as initial value in first iteration (no initializeData was called)
        BOOST_TEST(readData == 0);
      } else if (context.isNamed("SolverOne") && iterations == 0 && timewindow > 0) { // always use constant extrapolation in first iteration (from writeData of second participant at end previous window).
        BOOST_TEST(readData == readFunction(time));
      } else if (((context.isNamed("SolverOne") && iterations > 0) || context.isNamed("SolverTwo")) && timewindow == 0) {
        // first window is special, because of interpolation between zero initial data and data at end of first window.
        BOOST_TEST(readData == (readTime) / windowDt * readFunction(windowDt)); // self-made linear interpolation.
      } else if (((context.isNamed("SolverOne") && iterations > 0) ||           // SolverOne can only perform interpolation in later iterations
                  context.isNamed("SolverTwo"))                                 // Always use linear interpolation for SolverTwo
                 && timewindow > 0) {
        BOOST_TEST(readData == readFunction(readTime));
      } else {
        BOOST_TEST(false); // unreachable!
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    writeData = writeFunction(time);
    precice.writeScalarData(meshName, writeDataName, vertexID, writeData);
    precice.advance(maxDt);
    double maxDt = precice.getMaxTimeStepSize();
    currentDt    = dt > maxDt ? maxDt : dt;
    BOOST_CHECK(currentDt == windowDt); // no subcycling.
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
BOOST_AUTO_TEST_SUITE_END() // ParallelCoupling

#endif // PRECICE_NO_MPI
