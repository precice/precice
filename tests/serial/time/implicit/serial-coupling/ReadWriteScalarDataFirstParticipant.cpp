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
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataFirstParticipant)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  SolverInterface precice(context.name, context.config(), 0, 1);

  MeshID meshID;
  DataID writeDataID;
  DataID readDataID;

  // SolverOne prescribes these, thus SolverTwo expect these (we use "first-participant" as dt method)
  std::vector<std::vector<double>> timestepSizes{{1.0, 2.0, 1.0}, {2.0, 1.0, 2.0}, {3.0, 2.5, 3.0}};

  // max number of iterations in implicit coupling
  int maxIterations = 3;

  // some dummy values, to check the actual values is not the point of this test,
  // but more to test whether reading / writing is possible at all
  double expectedDataValue = 2.5;
  double actualDataValue   = -1.0;

  if (context.isNamed("SolverOne")) {
    auto meshID      = "MeshOne";
    auto writeDataID = "DataOne"; //  meshID
    auto readDataID  = "DataTwo"; //  meshID
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto meshID      = "MeshTwo";
    auto writeDataID = "DataTwo"; //  meshID
    auto readDataID  = "DataOne"; //  meshID
  }

  VertexID vertexID = precice.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
  double   dt       = precice.initialize();

  if (precice.requiresWritingCheckpoint()) {
    // do nothing
  }

  double startOfWindowTime = 0;
  double timeInWindow      = 0;
  double totalTime         = 6; // max-time from config

  for (auto iterationSizes : timestepSizes) {
    for (int it = 0; it < maxIterations; it++) {
      actualDataValue = -1; // reset value.
      BOOST_TEST(precice.isCouplingOngoing());
      precice.writeScalarData(meshID, writeDataID, vertexID, expectedDataValue);

      if (context.isNamed("SolverOne")) {
        dt = precice.advance(iterationSizes.at(it));
        timeInWindow += iterationSizes.at(it);
      } else if (context.isNamed("SolverTwo")) {
        BOOST_TEST(dt == iterationSizes.at(it));
        dt = precice.advance(dt);
      }

      if (precice.requiresReadingCheckpoint()) {
        timeInWindow = 0;
      }
      if (precice.isTimeWindowComplete()) {
        startOfWindowTime += timeInWindow;
        timeInWindow = 0;
      }

      if (context.isNamed("SolverOne")) {
        // Check remainder of simulation time
        BOOST_TEST(dt == totalTime - startOfWindowTime);
      }

      if (precice.requiresWritingCheckpoint()) {
        // do nothing
      }

      precice.readScalarData(meshID, readDataID, vertexID, actualDataValue);
      BOOST_TEST(actualDataValue == expectedDataValue);
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
