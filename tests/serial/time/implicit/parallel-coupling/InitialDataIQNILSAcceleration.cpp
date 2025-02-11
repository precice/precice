#ifndef PRECICE_NO_MPI

#include <math.h>
#include <precice/precice.hpp>
#include <vector>
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Simple test to ensure that underrelaxation is applied to every  timestep.
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(InitialDataIQNILSAcceleration)
{
  PRECICE_TEST();

  Participant precice(context.name, context.config(), 0, 1);

  MeshID meshID;
  DataID writeDataID;
  DataID readDataID;

  typedef double (*DataFunction)(double);

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) 2 + (2 * t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) 300 + (300 * t);
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
  double   v0[]      = {0, 0, 0};
  VertexID vertexID;
  vertexID = precice.setMeshVertex(meshName, v0);

  int    timestep           = 0;
  double time               = 0;
  int    timestepCheckpoint = timestep;
  double timeCheckpoint     = time;
  int    iterations         = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  double maxDt     = precice.getMaxTimeStepSize();
  double windowDt  = maxDt;
  double dt        = maxDt;
  double currentDt = dt; // Timestep length used by solver

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      timeCheckpoint     = time;
      timestepCheckpoint = timestep;
      iterations         = 0;
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    currentDt = dt > maxDt ? maxDt : dt;

    // only check data at beginning of window and ensure it does not change
    precice.readData(meshName, readDataName, {&vertexID, 1}, 0, {&readData, 1});
    BOOST_TEST(readData == readFunction(time));

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
  }

  precice.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
