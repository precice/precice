#include <boost/test/tools/context.hpp>
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
 * @brief Test to run a simple serial coupling where the first participant prescribes the time window size. Maximum simulation time is defined via max-time-windows, not max-time.
 *
 * Ensures that time window sizes are passed correctly and that reading and writing is possible.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataFirstParticipantFixedWindows)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  // SolverOne prescribes these, thus SolverTwo expect these (we use "first-participant" as dt method)
  std::vector<std::vector<double>> timestepSizes{{1.0, 2.0, 1.0}, {2.0, 1.0, 2.0}, {3.0, 2.5, 3.0}};

  // max number of iterations in implicit coupling
  int maxIterations = 3;

  // some dummy values, to check the actual values is not the point of this test,
  // but more to test whether reading / writing is possible at all
  const double expectedDataValue = 2.5;

  std::string meshName, writeDataName, readDataName;
  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "DataOne";
    readDataName  = "DataTwo";
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "DataTwo";
    readDataName  = "DataOne";
  }

  double   v0[]     = {0, 0, 0};
  VertexID vertexID = precice.setMeshVertex(meshName, v0);
  precice.initialize();
  double dt = precice.getMaxTimeStepSize();

  BOOST_TEST(precice.requiresWritingCheckpoint());

  int tw = 0;
  for (auto iterationSizes : timestepSizes) {
    for (int it = 0; it < maxIterations; it++) {
      BOOST_TEST_CONTEXT("tw " << tw << ", it " << it)
      {
        BOOST_TEST(precice.isCouplingOngoing());
        precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&expectedDataValue, 1});

        if (context.isNamed("SolverOne")) {
          dt = iterationSizes.at(it);
          precice.advance(dt);
          BOOST_TEST(precice.getMaxTimeStepSize() == std::numeric_limits<double>::max());
        } else if (context.isNamed("SolverTwo")) {
          BOOST_TEST(dt == iterationSizes.at(it));
          precice.advance(dt);
          dt = precice.getMaxTimeStepSize();
        }

        BOOST_TEST(precice.requiresReadingCheckpoint() == (it != 2));

        if (it == 2 && tw != 2) {
          BOOST_TEST(precice.requiresWritingCheckpoint());
        } else {
          BOOST_TEST(!precice.requiresWritingCheckpoint());
        }

        double receivedValue{-1};
        precice.readData(meshName, readDataName, {&vertexID, 1}, dt, {&receivedValue, 1});
        BOOST_TEST(receivedValue == expectedDataValue);
      }
    }
    ++tw;
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
