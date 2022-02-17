#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(PreciceTests)
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

  MeshID meshID;
  DataID writeDataID;
  DataID readDataID;

  typedef double (*DataFunction)(double, int);

  DataFunction dataOneFunction = [](double t, int idx) -> double {
    return (double) (2 + t + idx);
  };
  DataFunction dataTwoFunction = [](double t, int idx) -> double {
    return (double) (10 + t + idx);
  };
  DataFunction writeFunction;
  DataFunction readFunction;

  if (context.isNamed("SolverOne")) {
    meshID        = precice.getMeshID("MeshOne");
    writeDataID   = precice.getDataID("DataOne", meshID);
    writeFunction = dataOneFunction;
    readDataID    = precice.getDataID("DataTwo", meshID);
    readFunction  = dataTwoFunction;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshID        = precice.getMeshID("MeshTwo");
    writeDataID   = precice.getDataID("DataTwo", meshID);
    writeFunction = dataTwoFunction;
    readDataID    = precice.getDataID("DataOne", meshID);
    readFunction  = dataOneFunction;
  }

  int nVertices = 2;

  std::vector<VertexID> vertexIDs(nVertices, 0);
  std::vector<double>   writeData(nVertices, 0);
  std::vector<double>   readData(nVertices, 0);

  vertexIDs[0] = precice.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
  vertexIDs[1] = precice.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());

  int    nWindows   = 5; // perform 5 windows.
  double maxDt      = precice.initialize();
  double windowDt   = maxDt;
  int    timewindow = 0;
  int    timewindowCheckpoint;
  double dt        = maxDt; // Timestep length desired by solver
  double currentDt = dt;    // Timestep length used by solver
  double time      = timewindow * dt;
  double timeCheckpoint;
  double sampleDts[4] = {0.0, dt / 4.0, dt / 2.0, 3.0 * dt / 4.0};
  double readDts[4]   = {currentDt, currentDt, currentDt, currentDt};
  int    nSamples     = 4;
  int    iterations   = 0;
  double sampleDt; // dt relative to timestep start, where we are sampling
  double readDt;   // dt relative to timestep start for readTime
  double readTime; // time where we are reading from the reference solution
  if (precice.isActionRequired(precice::constants::actionWriteInitialData())) {
    for (int i = 0; i < nVertices; i++) {
      writeData[i] = writeFunction(time, i);
      precice.writeScalarData(writeDataID, vertexIDs[i], writeData[i]);
    }
    precice.markActionFulfilled(precice::constants::actionWriteInitialData());
  }

  precice.initializeData();

  while (precice.isCouplingOngoing()) {
    if (precice.isActionRequired(precice::constants::actionWriteIterationCheckpoint())) {
      timeCheckpoint       = time;
      timewindowCheckpoint = timewindow;
      iterations           = 0;
      precice.markActionFulfilled(precice::constants::actionWriteIterationCheckpoint());
    }
    BOOST_TEST(precice.isReadDataAvailable());
    BOOST_TEST(readData.size() == nVertices);
    for (int i = 0; i < nVertices; i++) {
      for (int j = 0; j < nSamples; j++) {
        sampleDt = sampleDts[j];
        readDt   = readDts[j];
        readTime = time + readDt;
        if (precice.isReadDataAvailable()) {
          precice.readScalarData(readDataID, vertexIDs[i], sampleDt, readData[i]);
        }
        if (iterations == 0) {
          BOOST_TEST(readData[i] == readFunction(time, i));
        } else if (iterations > 0) {
          BOOST_TEST(readData[i] == readFunction(readTime, i));
        } else {
          BOOST_TEST(false); // unreachable!
        }
      }
    }
    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    timewindow++;
    for (int i = 0; i < nVertices; i++) {
      writeData[i] = writeFunction(time, i);
    }

    if (precice.isWriteDataRequired(currentDt)) {
      BOOST_TEST(writeData.size() == nVertices);
      for (int i = 0; i < nVertices; i++) {
        writeData[i] = writeFunction(time, i);
        precice.writeScalarData(writeDataID, vertexIDs[i], writeData[i]);
      }
    }
    maxDt = precice.advance(currentDt);
    if (precice.isActionRequired(precice::constants::actionReadIterationCheckpoint())) {
      time       = timeCheckpoint;
      timewindow = timewindowCheckpoint;
      iterations++;
      precice.markActionFulfilled(precice::constants::actionReadIterationCheckpoint());
    }
    currentDt = dt > maxDt ? maxDt : dt;
  }

  precice.finalize();
  BOOST_TEST(timewindow == nWindows);
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
