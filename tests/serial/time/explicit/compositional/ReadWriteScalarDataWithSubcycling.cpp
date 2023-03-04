#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(Compositional)

/**
 * @brief Test to run a simple coupling with subcycling.
 *
 * Ensures that each time step provides its own data, but preCICE only exchanges data at the end of the window.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithSubcycling)
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
    writeDataName = "DataOne"; //  meshName
    writeFunction = dataOneFunction;
    readDataName  = "DataTwo"; //  meshName
    readFunction  = dataTwoFunction;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "DataTwo"; //  meshName
    writeFunction = dataTwoFunction;
    readDataName  = "DataOne"; //  meshName
    readFunction  = dataOneFunction;
  }

  double writeData, readData;

  VertexID vertexID = precice.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());

  int    nSubsteps  = 4; // perform subcycling on solvers. 4 steps happen in each window.
  int    nWindows   = 5; // perform 5 windows.
  int    timestep   = 0;
  double time       = 0;
  int    timewindow = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeScalarData(meshName, writeDataName, vertexID, writeData);
  }

  double maxDt = precice.initialize();
  BOOST_TEST(maxDt == 2.0); // use window size != 1.0 to be able to detect more possible bugs
  double windowDt      = maxDt;
  double dt            = windowDt / (nSubsteps - 0.5);                 // Solver tries to do a timestep of size 4/7
  double expectedDts[] = {4.0 / 7.0, 4.0 / 7.0, 4.0 / 7.0, 2.0 / 7.0}; // If solver uses timestep size of 4/7, fourth step will be restricted to 2/7 via preCICE steering to fit into the window.
  double currentDt     = dt > maxDt ? maxDt : dt;                      // determine actual timestep length; must fit into remaining time in window

  while (precice.isCouplingOngoing()) {
    double readTime = timewindow * windowDt; // both solvers lag one window behind for parallel-explicit coupling.

    precice.readScalarData(meshName, readDataName, vertexID, readData);
    BOOST_TEST(readData == readFunction(readTime));

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    BOOST_TEST(currentDt == expectedDts[timestep % nSubsteps]);
    time += currentDt;

    writeData = writeFunction(time);
    precice.writeScalarData(meshName, writeDataName, vertexID, writeData);

    maxDt     = precice.advance(currentDt);
    currentDt = dt > maxDt ? maxDt : dt;
    timestep++;
    if (precice.isTimeWindowComplete()) {
      timewindow++;
    }
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
