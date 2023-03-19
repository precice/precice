#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(ParallelCoupling)

/**
 * @brief Test to run a simple coupling with sampling from a zeroth order waveform; no subcycling.
 *
 * Provides a dt argument to the read function.
 */

BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSamplingZero)
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

  double   writeData, readData;
  VertexID vertexID = precice.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());

  int    nWindows   = 5; // perform 5 windows.
  int    timewindow = 0;
  int    timewindowCheckpoint;
  double time = 0;
  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeScalarData(meshName, writeDataName, vertexID, writeData);
  }

  double maxDt     = precice.initialize();
  double dt        = maxDt; // Timestep length desired by solver
  double currentDt = dt;    // Timestep length used by solver
  double timeCheckpoint;
  double sampleDts[4] = {0.0, dt / 4.0, dt / 2.0, 3.0 * dt / 4.0};
  double readDts[4]   = {0.0, currentDt, currentDt, currentDt};
  int    nSamples     = 4;
  int    iterations   = 0;
  double sampleDt; // dt relative to timestep start, where we are sampling
  double readDt;   // dt relative to timestep start for readTime
  double readTime; // time where we are reading from the reference solution

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      timeCheckpoint       = time;
      timewindowCheckpoint = timewindow;
      iterations           = 0;
    }

    for (int j = 0; j < nSamples; j++) {
      sampleDt = sampleDts[j];
      readDt   = readDts[j];
      readTime = time + readDt;

      precice.readScalarData(meshName, readDataName, vertexID, sampleDt, readData);

      if (iterations == 0) {
        BOOST_TEST(readData == readFunction(time));
      } else if (iterations > 0) {
        BOOST_TEST(readData == readFunction(readTime));
      } else {
        BOOST_TEST(false); // unreachable!
      }
    }
    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    timewindow++;
    writeData = writeFunction(time);
    precice.writeScalarData(meshName, writeDataName, vertexID, writeData);
    maxDt = precice.advance(currentDt);
    if (precice.requiresReadingCheckpoint()) {
      time       = timeCheckpoint;
      timewindow = timewindowCheckpoint;
      iterations++;
    }
    currentDt = dt > maxDt ? maxDt : dt;
  }

  precice.finalize();
  BOOST_TEST(timewindow == nWindows);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // ParallelCoupling

#endif // PRECICE_NO_MPI
