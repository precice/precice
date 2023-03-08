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
 * @brief Test to run a simple serial coupling where the first participant prescribes the time window size.
 *
 * Ensures that time window sizes are passed correctly and that reading and writing is possible.
 *
 * Tests how Waveform relaxation and changing time window sizes (inside one window) interact. See https://github.com/precice/precice/issues/1570#issuecomment-1436067426.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataFirstParticipantChangingDt)
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

  // SolverOne prescribes these, thus SolverTwo expect these (we use "first-participant" as dt method)
  std::vector<std::vector<double>> timestepSizes{{1.0, 2.0, 1.0}, {2.0, 1.0, 2.0}, {3.0, 2.5, 3.0}};

  // max number of iterations in implicit coupling
  int maxIterations = 3;

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

  double writeData = 0;
  double readData  = 0;

  VertexID vertexID = precice.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());
  if (precice.requiresInitialData()) {
    precice.writeScalarData(meshName, writeDataName, vertexID, writeFunction(0));
  }
  double preciceDt = precice.initialize();
  double solverDt  = 0.5;
  double dt;
  double startOfWindowTime = 0;
  double timeInWindow      = 0;
  double totalTime         = 6; // max-time from config
  double expectedDataValue, actualDataValue;
  int    it     = 0;
  int    window = 0;

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      // do nothing
    }

    actualDataValue = -1; // reset value.

    if (context.isNamed("SolverOne")) {
      solverDt = timestepSizes.at(window).at(it); // SolverOne uses varying dt
    } else if (context.isNamed("SolverTwo")) {
      solverDt = 0.5; // SolverTwo uses fixed dt
    }

    // explicitly read at end of step
    if (context.isNamed("SolverOne")) {
      // @todo window end is not well defined for SolverOne! See https://github.com/precice/precice/issues/1570#issuecomment-1443063091
    } else {
      precice.readScalarData(meshName, readDataName, vertexID, solverDt, actualDataValue);
    }

    dt = std::min({solverDt, preciceDt});
    timeInWindow += dt;

    if (context.isNamed("SolverOne")) {
      // @todo Not testing, because we actually receive the value at the end of the window from the last iteration of solver two. This is not necessary consistent with the end of this iteration of solver one. See https://github.com/precice/precice/issues/1570#issuecomment-1443063091
    } else {
      expectedDataValue = readFunction(startOfWindowTime + timeInWindow);
      BOOST_TEST(actualDataValue == expectedDataValue);
    }
    precice.writeScalarData(meshName, writeDataName, vertexID, writeFunction(startOfWindowTime + timeInWindow));
    preciceDt = precice.advance(dt);

    if (precice.requiresReadingCheckpoint()) {
      timeInWindow = 0;
      it++;
    }

    if (precice.isTimeWindowComplete()) {
      startOfWindowTime += timeInWindow;
      timeInWindow = 0;
      it           = 0;
      window++;
    }

    if (context.isNamed("SolverOne")) {
      BOOST_TEST(preciceDt == totalTime - startOfWindowTime);
    } else if (context.isNamed("SolverTwo") && precice.isCouplingOngoing()) {
      BOOST_TEST(preciceDt == timestepSizes.at(window).at(it) - timeInWindow);
    }
  }

  BOOST_TEST(not precice.isCouplingOngoing());
  precice.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
