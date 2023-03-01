#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(TestImplicit)
{
  /// Test simple coupled simulation with coupling iterations.
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  double state              = 0.0;
  double checkpoint         = 0.0;
  int    iterationCount     = 0;
  double initialStateChange = 5.0;
  double stateChange        = initialStateChange;
  int    computedTimesteps  = 0;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("SolverOne")) {
    auto   meshID = "Square";
    double pos[3];
    // Set mesh positions
    pos[0] = 0.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
    interface.setMeshVertex(meshID, pos);
    pos[0] = 1.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
    interface.setMeshVertex(meshID, pos);
    pos[0] = 1.0;
    pos[1] = 1.0;
    pos[2] = 0.0;
    interface.setMeshVertex(meshID, pos);
    pos[0] = 0.0;
    pos[1] = 1.0;
    pos[2] = 0.0;
    interface.setMeshVertex(meshID, pos);

    double maxDt = interface.initialize();
    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        checkpoint     = state;
        iterationCount = 1;
      }
      if (interface.requiresReadingCheckpoint()) {
        state = checkpoint;
      }
      iterationCount++;
      stateChange = initialStateChange / (double) iterationCount;
      state += stateChange;
      maxDt = interface.advance(maxDt);
      if (interface.isTimeWindowComplete()) {
        computedTimesteps++;
      }
    }
    interface.finalize();
    BOOST_TEST(computedTimesteps == 4);
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto   meshID = "SquareTwo";
    double pos[3];
    // Set mesh positions
    pos[0] = 0.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
    interface.setMeshVertex(meshID, pos);
    double maxDt = interface.initialize();
    while (interface.isCouplingOngoing()) {
      if (interface.requiresWritingCheckpoint()) {
        checkpoint     = state;
        iterationCount = 1;
      }
      if (interface.requiresReadingCheckpoint()) {
        state = checkpoint;
        iterationCount++;
      }
      stateChange = initialStateChange / (double) iterationCount;
      state += stateChange;
      maxDt = interface.advance(maxDt);
      if (interface.isTimeWindowComplete()) {
        computedTimesteps++;
      }
    }
    interface.finalize();
    BOOST_TEST(computedTimesteps == 4);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
