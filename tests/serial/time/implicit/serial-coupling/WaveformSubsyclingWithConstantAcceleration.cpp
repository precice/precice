#ifndef PRECICE_NO_MPI

#include <math.h>
#include <precice/SolverInterface.hpp>
#include <vector>
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Simple test to ensure that acceleration is applied to every  timestep.
 */
BOOST_AUTO_TEST_CASE(WaveformSubsyclingWithConstantAcceleration)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  SolverInterface precice(context.name, context.config(), 0, 1);

  MeshID meshID;
  DataID writeDataID;
  DataID readDataID;

  typedef double (*DataFunction)(double);

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 * t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (3 * t);
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

  double   writeData = 0;
  double   readData  = 0;
  VertexID vertexID  = precice.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());

  int    nSubsteps          = 7; // perform subcycling on solvers. 4 steps happen in each window.
  int    nWindows           = 5; // perform 5 windows.
  int    timestep           = 0;
  double time               = 0;
  int    timestepCheckpoint = timestep;
  double timeCheckpoint     = time;
  int    iterations         = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeScalarData(writeDataID, vertexID, writeData);
  }

  double maxDt     = precice.initialize();
  double windowDt  = maxDt;
  double dt        = windowDt / nSubsteps; // Timestep length desired by solver. E.g. 4 steps  with size 1/4
  double currentDt = dt;                   // Timestep length used by solver

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      timeCheckpoint     = time;
      timestepCheckpoint = timestep;
      iterations         = 0;
    }

    double factor = 1 - pow((1 - 0.1), iterations);
    precice.readScalarData(readDataID, vertexID, maxDt, readData);
    if (context.isNamed("SolverOne")) { // in the first iteration of each window, use data from previous window.
      if (iterations == 0) {
        BOOST_TEST(readData == readFunction(timeCheckpoint));
      } else {
        BOOST_TEST(readData == factor * readFunction(time + maxDt));
      }
    } else { // in the following iterations, use data at the end of window.
      BOOST_TEST(readData == readFunction(time + maxDt));
    }

    precice.readScalarData(readDataID, vertexID, currentDt, readData);
    if (context.isNamed("SolverOne")) { // in the first iteration of each window, use data from previous window.
      if (iterations == 0) {
        BOOST_TEST(readData == readFunction(timeCheckpoint));
      } else {
        BOOST_TEST(readData == factor * readFunction(time + currentDt));
      }
    } else { // in the following iterations, use data at the end of window.
      BOOST_TEST(readData == readFunction(time + currentDt));
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    timestep++;
    writeData = writeFunction(time);
    precice.writeScalarData(writeDataID, vertexID, writeData);
    maxDt = precice.advance(currentDt);
    if (precice.requiresReadingCheckpoint()) {
      time     = timeCheckpoint;
      timestep = timestepCheckpoint;
      iterations++;
    }
    currentDt = dt > maxDt ? maxDt : dt;
  }

  precice.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // ParallelCoupling

#endif // PRECICE_NO_MPI
