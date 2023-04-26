#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
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

  SolverInterface precice(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    precice.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    precice.setMeshVertex(meshName, Eigen::Vector3d(1.0, 0.0, 0.0).data());
    double maxDt     = precice.initialize();
    int    timestep  = 0;
    double dt        = maxDt / 2.0; // Time step size desired by solver
    double currentDt = dt;          // Time step size used by solver
    while (precice.isCouplingOngoing()) {
      maxDt     = precice.advance(currentDt);
      currentDt = dt > maxDt ? maxDt : dt;
      timestep++;
    }
    precice.finalize();
    BOOST_TEST(timestep == 20);
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto meshName = "Test-Square";
    precice.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    precice.setMeshVertex(meshName, Eigen::Vector3d(1.0, 0.0, 0.0).data());
    double maxDt     = precice.initialize();
    int    timestep  = 0;
    double dt        = maxDt / 3.0; // Time step size desired by solver
    double currentDt = dt;          // Time step size used by solver
    while (precice.isCouplingOngoing()) {
      maxDt     = precice.advance(currentDt);
      currentDt = dt > maxDt ? maxDt : dt;
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
