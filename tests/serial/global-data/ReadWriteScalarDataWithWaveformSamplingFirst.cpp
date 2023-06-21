#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(GlobalData)
/**
 * @brief Test to run a simple coupling with first order waveform subcycling
 * also applied to a global data.
 *
 * Provides a dt argument to the read function. A first order waveform is used.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSamplingFirst)
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
  DataFunction globalDataOneFunction = [](double t) -> double {
    return (double) (12 + t);
  };
  DataFunction globalDataTwoFunction = [](double t) -> double {
    return (double) (20 + t);
  };
  DataFunction writeFunction;
  DataFunction readFunction;
  DataFunction writeGlobalFunction;
  DataFunction readGlobalFunction;

  std::string meshName, writeDataName, readDataName, writeGlobalDataName, readGlobalDataName;
  if (context.isNamed("SolverOne")) {
    meshName            = "MeshOne";
    writeDataName       = "DataOne";
    writeFunction       = dataOneFunction;
    readDataName        = "DataTwo";
    readFunction        = dataTwoFunction;
    writeGlobalDataName = "GlobalDataOne";
    writeGlobalFunction = globalDataOneFunction;
    readGlobalDataName  = "GlobalDataTwo";
    readGlobalFunction  = globalDataTwoFunction;

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName            = "MeshTwo";
    writeDataName       = "DataTwo";
    writeFunction       = dataTwoFunction;
    readDataName        = "DataOne";
    readFunction        = dataOneFunction;
    writeGlobalDataName = "GlobalDataTwo";
    writeGlobalFunction = globalDataTwoFunction;
    readGlobalDataName  = "GlobalDataOne";
    readGlobalFunction  = globalDataOneFunction;
  }

  double   writeData       = 0;
  double   readData        = 0;
  double   writeGlobalData = 0;
  double   readGlobalData  = 0;
  double   v0[]            = {0, 0, 0};
  VertexID vertexID        = precice.setMeshVertex(meshName, v0);

  int    nWindows        = 5; // perform 5 windows.
  int    timestep        = 0;
  int    timewindow      = 0;
  double windowStartTime = 0;
  int    windowStartStep = 0;
  int    nSamples        = 4;
  int    iterations      = 0;
  double time            = 0;
  double readTime; // time where we are reading
  double sampleDt; // dt relative to timestep start, where we are sampling

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    writeGlobalData = writeGlobalFunction(time);
    precice.writeGlobalData(writeGlobalDataName, {&writeGlobalData, 1});
  }

  precice.initialize();
  double maxDt        = precice.getMaxTimeStepSize();
  double windowDt     = maxDt;
  double dt           = maxDt; // time step size desired by solver
  double currentDt    = dt;    // time step size used by solver
  double sampleDts[4] = {0.0, dt / 4.0, dt / 2.0, 3.0 * dt / 4.0};

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      windowStartTime = time;
      windowStartStep = timestep;
    }
    for (int j = 0; j < nSamples; j++) {
      sampleDt = sampleDts[j];
      readTime = time + sampleDt;
      precice.readData(meshName, readDataName, {&vertexID, 1}, sampleDt, {&readData, 1});
      precice.readGlobalData(readGlobalDataName, sampleDt, {&readGlobalData, 1});

      if (context.isNamed("SolverOne") && iterations == 0) { // first participant always uses constant extrapolation in first iteration (from initializeData or writeData of second participant at end previous window).
        BOOST_TEST(readData == readFunction(time));
        BOOST_TEST(readGlobalData == readGlobalFunction(time));
      } else { // second participant always uses linear interpolation in later windows (additionally available writeData of first participant at end of this window).
        BOOST_TEST(readData == readFunction(readTime));
        BOOST_TEST(readGlobalData == readGlobalFunction(readTime));
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    writeData       = writeFunction(time);
    writeGlobalData = writeGlobalFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.writeGlobalData(writeGlobalDataName, {&writeGlobalData, 1});
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
BOOST_AUTO_TEST_SUITE_END() // GlobalData

#endif // PRECICE_NO_MPI
