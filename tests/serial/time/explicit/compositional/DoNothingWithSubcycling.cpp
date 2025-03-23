#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(Compositional)

/**
 * @brief Test to run a "do nothing" compositional coupling with subcycling solvers.
 *
 */
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank))
BOOST_AUTO_TEST_CASE(DoNothingWithSubcycling)
{
  PRECICE_TEST();

  Participant precice(context.name, context.config(), 0, 1);

  int nSubsteps; // let three solvers use different time step sizes

  std::string meshName, writeDataName;
  if (context.isNamed("SolverOne")) {
    meshName  = "MeshOne";
    nSubsteps = 1;
  } else if (context.isNamed("SolverTwo")) {
    meshName  = "MeshTwo";
    nSubsteps = 2;
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    meshName  = "MeshThree";
    nSubsteps = 3;
  }

  double v0[] = {0, 0, 0};
  precice.setMeshVertex(meshName, v0);

  int totalSolves             = 0;
  int totalCompletedTimesteps = 0;
  int timestepsInThisWindow   = 0;

  precice.initialize();
  double maxDt     = precice.getMaxTimeStepSize();
  double windowDt  = maxDt;
  double dt        = windowDt / nSubsteps; // time step size desired by solver. E.g. 2 steps with size 1/2
  double currentDt = dt;                   // time step size used by solver

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      totalCompletedTimesteps += timestepsInThisWindow;
      timestepsInThisWindow = 0;
    }
    precice.advance(currentDt);
    double maxDt = precice.getMaxTimeStepSize();
    currentDt    = dt > maxDt ? maxDt : dt;
    totalSolves++;
    timestepsInThisWindow++;
    if (precice.requiresReadingCheckpoint()) {
      timestepsInThisWindow = 0;
    }
  }

  totalCompletedTimesteps += timestepsInThisWindow;
  precice.finalize();

  if (context.isNamed("SolverOne")) {
    BOOST_TEST(totalCompletedTimesteps == 5);
    BOOST_TEST(totalSolves == 5);
  } else if (context.isNamed("SolverTwo")) {
    BOOST_TEST(totalCompletedTimesteps == 10);
    BOOST_TEST(totalSolves == 10);
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    BOOST_TEST(totalCompletedTimesteps == 15);
    BOOST_TEST(totalSolves == 15);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // Compositional

#endif // PRECICE_NO_MPI
