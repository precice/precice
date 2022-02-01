#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to run a simple coupling with first order waveform subcycling.
 * 
 * Provides a dt argument to the read function. A first order waveform is used.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingFirst)
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
  double dt        = windowDt / nSubsteps; // Timestep length desired by solver. E.g. 4 steps  with size 1/4
  dt += windowDt / nSubsteps / nSubsteps;  // increase timestep such that we get a non-matching subcycling. E.g. 3 step with size 5/16 and 1 step with size 1/16.
  double currentDt = dt;                   // Timestep length used by solver
  double time      = timestep * dt;

  if (precice.isActionRequired(precice::constants::actionWriteInitialData())) {
    for (int i = 0; i < nVertices; i++) {
      writeData[i] = writeFunction(time, i);
      precice.writeScalarData(writeDataID, vertexIDs[i], writeData[i]);
    }
    precice.markActionFulfilled(precice::constants::actionWriteInitialData());
  }

  precice.initializeData();

  while (precice.isCouplingOngoing()) {
    double readTime;
    if (context.isNamed("SolverOne")) {
      readTime = time - windowDt + currentDt; // solver one lags one window behind solver two.
    } else {
      readTime = time + currentDt;
    }

    BOOST_TEST(readData.size() == nVertices);
    BOOST_TEST(precice.isReadDataAvailable());
    for (int i = 0; i < nVertices; i++) {
      if (precice.isReadDataAvailable()) {
        precice.readScalarData(readDataID, vertexIDs[i], currentDt, readData[i]);
      }
      if (timestep < nSubsteps) { // in the first window, we only have one sample of data. Therefore constant interpolation
        if (context.isNamed("SolverOne")) {
          BOOST_TEST(readData[i] == readFunction(0, i));
        } else {
          BOOST_TEST(readData[i] == readFunction(windowDt, i));
        }
      } else { // in the following windows we have two samples of data. Therefore linear interpolation
        BOOST_TEST(readData[i] == readFunction(readTime, i));
      }
      if (precice.isReadDataAvailable()) {
        precice.readScalarData(readDataID, vertexIDs[i], currentDt / 2, readData[i]);
      }
      if (timestep < nSubsteps) { // in the first window, we only have one sample of data. Therefore constant interpolation
        if (context.isNamed("SolverOne")) {
          BOOST_TEST(readData[i] == readFunction(0, i));
        } else {
          BOOST_TEST(readData[i] == readFunction(windowDt, i));
        }
      } else { // in the following windows we have two samples of data. Therefore linear interpolation
        BOOST_TEST(readData[i] == readFunction(readTime - currentDt / 2, i));
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
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
    maxDt     = precice.advance(currentDt);
    currentDt = dt > maxDt ? maxDt : dt;
    timestep++;
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
