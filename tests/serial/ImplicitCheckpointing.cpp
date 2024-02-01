#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(ImplicitCheckpointing)
{
  /// Test simple implicit coupling with checkpointing. Checks correct tracking of time, see https://github.com/precice/precice/pull/1704.
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  double state              = 0.0;
  double checkpoint         = 0.0;
  int    iterationCount     = 0;
  double initialStateChange = 5.0;
  double stateChange        = initialStateChange;
  int    computedTimesteps  = 0;

  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("SolverOne")) {
    auto   meshName = "Square";
    double pos[3];
    // Set mesh positions
    pos[0] = 0.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
    interface.setMeshVertex(meshName, pos);
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto   meshName = "SquareTwo";
    double pos[3];
    // Set mesh positions
    pos[0] = 0.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
    interface.setMeshVertex(meshName, pos);
  }

  interface.initialize();
  double maxDt = interface.getMaxTimeStepSize();

  BOOST_TEST(interface.isCouplingOngoing());
  BOOST_TEST(interface.requiresWritingCheckpoint());
  BOOST_TEST(not interface.requiresReadingCheckpoint());

  interface.advance(maxDt); // finish first iteration of first window

  BOOST_TEST(interface.isCouplingOngoing());
  BOOST_TEST(not interface.requiresWritingCheckpoint());
  BOOST_TEST(interface.requiresReadingCheckpoint());

  interface.advance(maxDt); // finish second iteration of first window

  BOOST_TEST(interface.isCouplingOngoing());
  BOOST_TEST(not interface.requiresWritingCheckpoint());
  BOOST_TEST(interface.requiresReadingCheckpoint());

  interface.advance(maxDt); // finish third and last iteration of first window

  BOOST_TEST(interface.isCouplingOngoing());
  BOOST_TEST(interface.requiresWritingCheckpoint());
  BOOST_TEST(not interface.requiresReadingCheckpoint());

  interface.advance(maxDt); // finish first iteration of second window

  BOOST_TEST(interface.isCouplingOngoing());
  BOOST_TEST(not interface.requiresWritingCheckpoint());
  BOOST_TEST(interface.requiresReadingCheckpoint());

  interface.advance(maxDt); // finish second iteration of second window

  BOOST_TEST(interface.isCouplingOngoing());
  BOOST_TEST(not interface.requiresWritingCheckpoint());
  BOOST_TEST(interface.requiresReadingCheckpoint());

  interface.advance(maxDt); // finish third and last iteration of second window

  BOOST_TEST(not interface.isCouplingOngoing());
  BOOST_TEST(not interface.requiresWritingCheckpoint());
  BOOST_TEST(not interface.requiresReadingCheckpoint());

  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
