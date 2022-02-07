#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to run a simple coupling with zeroth order waveform subcycling.
 * 
 * Provides a dt argument to the read function, but since a zeroth order waveform is used the result should be identical to the case without waveform relaxation
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingZero)
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

  int nVertices = 1;

  std::vector<VertexID> vertexIDs(nVertices, 0);
  std::vector<double>   writeData(nVertices, 0);
  std::vector<double>   readData(nVertices, 0);

  vertexIDs[0] = precice.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());

  int    nSubsteps = 4; // perform subcycling on solvers. 4 steps happen in each window.
  int    nWindows  = 5; // perform 5 windows.
  double maxDt     = precice.initialize();
  double windowDt  = maxDt;
  int    timestep  = 0;
  int    timestepCheckpoint;
  double dt = windowDt / nSubsteps;       // Timestep length desired by solver. E.g. 4 steps  with size 1/4
  dt += windowDt / nSubsteps / nSubsteps; // increase timestep such that we get a non-matching subcycling. E.g. 3 step with size 5/16 and 1 step with size 1/16.
  double currentDt = dt;                  // Timestep length used by solver
  double time      = timestep * dt;
  double timeCheckpoint;
  int    iterations;

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
      timeCheckpoint     = time;
      timestepCheckpoint = timestep;
      iterations         = 0;
      precice.markActionFulfilled(precice::constants::actionWriteIterationCheckpoint());
    }
    double readTime;
    readTime = timeCheckpoint + windowDt;

    BOOST_TEST(readData.size() == nVertices);
    BOOST_TEST(precice.isReadDataAvailable());
    for (int i = 0; i < nVertices; i++) {
      if (precice.isReadDataAvailable()) {
        precice.readScalarData(readDataID, vertexIDs[i], currentDt, readData[i]);
      }
      if (context.isNamed("SolverOne") && iterations == 0) { // in the first iteration of each window, use data from previous window.
        BOOST_TEST(readData[i] == readFunction(timeCheckpoint, i));
      } else { // in the following iterations, use data at the end of window.
        BOOST_TEST(readData[i] == readFunction(readTime, i));
      }
      if (precice.isReadDataAvailable()) {
        precice.readScalarData(readDataID, vertexIDs[i], currentDt / 2, readData[i]);
      }
      if (context.isNamed("SolverOne") && iterations == 0) { // in the first iteration of each window, use data from previous window.
        BOOST_TEST(readData[i] == readFunction(timeCheckpoint, i));
      } else { // in the following iterations, use data at the end of window.
        BOOST_TEST(readData[i] == readFunction(readTime, i));
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    timestep++;
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
      time     = timeCheckpoint;
      timestep = timestepCheckpoint;
      iterations++;
      precice.markActionFulfilled(precice::constants::actionReadIterationCheckpoint());
    }
    currentDt = dt > maxDt ? maxDt : dt;
  }

  precice.finalize();
  BOOST_TEST(timestep == nWindows * nSubsteps);
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
