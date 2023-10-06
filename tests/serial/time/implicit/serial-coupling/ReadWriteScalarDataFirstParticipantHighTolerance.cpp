#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to run a simple serial coupling with a high min-timestep.
 *
 * Ensures that the feature min time step works correctly.
 *
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataFirstParticipantHighTolerance)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  typedef double (*DataFunction)(double);

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (10 + t * t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (2 + t * t);
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

  double writeData = 0;
  double readData  = 0;

  double   v0[]     = {0, 0, 0};
  VertexID vertexID = precice.setMeshVertex(meshName, v0);
  if (precice.requiresInitialData()) {
    double data[] = {writeFunction(0)};
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, data);
  }
  precice.initialize();
  double preciceDt = precice.getMaxTimeStepSize();
  double dt;
  double solverDt;
  double startOfWindowTime = 0;
  double timeInWindow      = 0;
  double totalTime         = 6; // max-time from config
  double expectedDataValue, actualDataValue;
  int    iterations     = 0;
  int    window         = 0;
  int    nbrTimeWindows = 0;

  if (context.isNamed("SolverOne")) {
    solverDt = 0.12; // non matching fixed dt that is bigger than min-timestep
  } else if (context.isNamed("SolverTwo")) {
    solverDt = 0.11; // non matching fixed dt that is bigger than min-timestep
  }

  while (precice.isCouplingOngoing()) {

    if (precice.requiresWritingCheckpoint()) {
      startOfWindowTime += timeInWindow;
      timeInWindow = 0;
      iterations   = 0;
      nbrTimeWindows++;
    }

    dt = solverDt > preciceDt ? preciceDt : solverDt;

    //It is enough to only check this in the first time step and in the second iteration
    if (iterations > 0 && timeInWindow == 0) {
      if (context.isNamed("SolverOne")) {
        precice.readData(meshName, readDataName, {&vertexID, 1}, preciceDt, {&actualDataValue, 1});
        BOOST_TEST(actualDataValue == readFunction(0.99 * nbrTimeWindows)); // Check that the end of the time window has the correct value

        precice.readData(meshName, readDataName, {&vertexID, 1}, 0.99, {&actualDataValue, 1});
        BOOST_TEST(actualDataValue == readFunction(0.99 * nbrTimeWindows)); // Check that the last sample from the solver has the correct value

        precice.readData(meshName, readDataName, {&vertexID, 1}, 0.33, {&actualDataValue, 1});
        BOOST_TEST(math::equals(actualDataValue, readFunction(0.33 + 0.99 * (nbrTimeWindows - 1)), 1e-6)); // We loose accuracy when sampling

        precice.readData(meshName, readDataName, {&vertexID, 1}, 0.95, {&actualDataValue, 1});
        BOOST_TEST(!math::equals(actualDataValue, readFunction(0.95 + 0.99 * (nbrTimeWindows - 1)), 1e-2)); // The waveform has changed here due to the added time step

      } else {

        precice.readData(meshName, readDataName, {&vertexID, 1}, preciceDt, {&actualDataValue, 1});
        BOOST_TEST(actualDataValue == readFunction(0.96 * nbrTimeWindows)); // Check that the end of the time window has the correct value

        precice.readData(meshName, readDataName, {&vertexID, 1}, 0.96, {&actualDataValue, 1});
        BOOST_TEST(actualDataValue == readFunction(0.96 * nbrTimeWindows)); // Check that the last sample from the solver has the correct value
      }
    }

    timeInWindow += dt;
    double data[] = {writeFunction(startOfWindowTime + timeInWindow)};
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, data);
    precice.advance(dt);

    preciceDt = precice.getMaxTimeStepSize();

    if (precice.requiresReadingCheckpoint()) {
      timeInWindow = 0;
      iterations++;
    }

    if (precice.isTimeWindowComplete()) {
      if (context.isNamed("SolverOne")) {
        BOOST_TEST(startOfWindowTime + timeInWindow == 0.96 * nbrTimeWindows); // Check that both solvers have stopped at the correct timestep
      } else {
        BOOST_TEST(startOfWindowTime + timeInWindow == 0.99 * nbrTimeWindows); // Check that both solvers have stopped at the correct timestep
      }
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
