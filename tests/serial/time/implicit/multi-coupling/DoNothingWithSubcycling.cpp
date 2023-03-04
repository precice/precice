#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(MultiCoupling)

/**
 * @brief Test to run a "do nothing" multi coupling with subcycling solvers.
 *
 */
BOOST_AUTO_TEST_CASE(DoNothingWithSubcycling)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));

  SolverInterface precice(context.name, context.config(), 0, 1);

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

  VertexID vertexID = precice.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());

  int totalSolves             = 0;
  int totalCompletedTimesteps = 0;
  int timestepsInThisWindow   = 0;

  double maxDt     = precice.initialize();
  double windowDt  = maxDt;
  double dt        = windowDt / nSubsteps; // Timestep length desired by solver. E.g. 2 steps with size 1/2
  double currentDt = dt;                   // Timestep length used by solver

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      totalCompletedTimesteps += timestepsInThisWindow;
      timestepsInThisWindow = 0;
    }
    maxDt     = precice.advance(currentDt);
    currentDt = dt > maxDt ? maxDt : dt;
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
    BOOST_TEST(totalSolves == 15);
  } else if (context.isNamed("SolverTwo")) {
    BOOST_TEST(totalCompletedTimesteps == 10);
    BOOST_TEST(totalSolves == 30);
  } else {
    BOOST_TEST(context.isNamed("SolverThree"));
    BOOST_TEST(totalCompletedTimesteps == 15);
    BOOST_TEST(totalSolves == 45);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // MultiCoupling

#endif // PRECICE_NO_MPI
