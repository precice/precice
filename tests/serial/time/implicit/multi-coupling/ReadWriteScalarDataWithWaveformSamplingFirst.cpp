#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(MultiCoupling)
/**
 * @brief Test to run a multi coupling with first order time interpolation.
 *
 * Provides a dt argument to the read function. A first order waveform is used.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSamplingFirst)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));

  SolverInterface precice(context.name, context.config(), 0, 1);

  std::string meshID, writeDataID;

  typedef double (*DataFunction)(double);
  std::vector<std::pair<std::string, DataFunction>> readDataPairs;

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 + t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (10 + t);
  };
  DataFunction dataThreeFunction = [](double t) -> double {
    return (double) (300 + t);
  };
  DataFunction writeFunction;

  if (context.isNamed("SolverOne")) {
    meshID         = "MeshOne";
    writeDataID    = "DataOne"; //  meshID
    writeFunction  = dataOneFunction;
    auto dataTwoId = "DataTwo"; //  meshID
    readDataPairs.push_back(std::make_pair(dataTwoId, dataTwoFunction));
    auto dataThreeId = "DataThree"; //  meshID
    readDataPairs.push_back(std::make_pair(dataThreeId, dataThreeFunction));
  } else if (context.isNamed("SolverTwo")) {
    meshID         = "MeshTwo";
    writeDataID    = "DataTwo"; //  meshID
    writeFunction  = dataTwoFunction;
    auto dataOneId = "DataOne"; //  meshID
    readDataPairs.push_back(std::make_pair(dataOneId, dataOneFunction));
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    meshID         = "MeshThree";
    writeDataID    = "DataThree"; //  meshID
    writeFunction  = dataThreeFunction;
    auto dataOneId = "DataOne"; //  meshID
    readDataPairs.push_back(std::make_pair(dataOneId, dataOneFunction));
  }

  double   writeData = 0;
  double   readData  = 0;
  VertexID vertexID  = precice.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());

  int    nWindows        = 5; // perform 5 windows.
  int    timestep        = 0;
  int    timewindow      = 0;
  double windowStartTime = 0;
  int    windowStartStep = 0;
  int    nSamples        = 4;
  int    iterations      = 0;
  double time            = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeScalarData(meshID, writeDataID, vertexID, writeData);
  }

  double maxDt        = precice.initialize();
  double windowDt     = maxDt;
  double dt           = maxDt; // Timestep length desired by solver
  double currentDt    = dt;    // Timestep length used by solver
  double sampleDts[4] = {0.0, dt / 4.0, dt / 2.0, 3.0 * dt / 4.0};
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

      for (auto &readDataPair : readDataPairs) {
        auto readDataID   = readDataPair.first;
        auto readFunction = readDataPair.second;

        precice.readScalarData(meshID, readDataID, vertexID, sampleDt, readData);

        if (iterations == 0) { // always use constant extrapolation in first iteration (from initialize or writeData of second participant at end previous window).
          BOOST_TEST(readData == readFunction(time));
        } else if (iterations > 0) { // use linear interpolation in later iterations (additionally available writeData of second participant at end of this window).
          BOOST_TEST(readData == readFunction(readTime));
        } else {
          BOOST_TEST(false); // unreachable!
        }
      }
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    writeData = writeFunction(time);
    precice.writeScalarData(meshID, writeDataID, vertexID, writeData);
    maxDt     = precice.advance(currentDt);
    currentDt = dt > maxDt ? maxDt : dt;
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
BOOST_AUTO_TEST_SUITE_END() // MultiCoupling

#endif // PRECICE_NO_MPI
