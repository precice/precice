#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to run a simple "do nothing" coupling with subcycling solvers.
 *
 */
BOOST_AUTO_TEST_CASE(DoNothingWithSubcycling)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  double v0[] = {0, 0, 0};
  double v1[] = {1, 0, 0};

  Participant precice(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    precice.setMeshVertex(meshName, v0);
    precice.setMeshVertex(meshName, v1);
    precice.initialize();
    double maxDt     = precice.getMaxTimeStepSize();
    int    timestep  = 0;
    double dt        = maxDt / 2.0; // Time step size desired by solver
    double currentDt = dt;          // Time step size used by solver
    while (precice.isCouplingOngoing()) {
      precice.advance(currentDt);
      maxDt     = precice.getMaxTimeStepSize();
      currentDt = dt > maxDt ? maxDt : dt;
      timestep++;
    }
    precice.finalize();
    BOOST_TEST(timestep == 20);
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto meshName = "Test-Square";
    precice.setMeshVertex(meshName, v0);
    precice.setMeshVertex(meshName, v1);
    precice.initialize();
    double maxDt     = precice.getMaxTimeStepSize();
    int    timestep  = 0;
    double dt        = maxDt / 3.0; // Time step size desired by solver
    double currentDt = dt;          // Time step size used by solver
    while (precice.isCouplingOngoing()) {
      maxDt     = precice.getMaxTimeStepSize();
      currentDt = dt > maxDt ? maxDt : dt;
      precice.advance(currentDt);
      timestep++;
    }
    precice.finalize();
    BOOST_TEST(timestep == 30);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
