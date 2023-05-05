#include <boost/test/unit_test_log.hpp>
#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "precice/SolverInterface.hpp"

/// Test to run simple "do nothing" coupling between two solvers.
void runTestExplicit(std::string const &configurationFileName, TestContext const &context)
{
  BOOST_TEST_MESSAGE("Config: " << configurationFileName);

  int    timesteps = 0;
  double time      = 0.0;

  SolverInterface couplingInterface(context.name, configurationFileName, 0, 1);

  double pos[] = {0, 0, 0, 1, 1, 1};
  int    vids[2];

  // was necessary to replace pre-defined geometries
  if (context.isNamed("SolverOne")) {
    auto meshName = "MeshOne";
    BOOST_REQUIRE(couplingInterface.hasMesh(meshName));
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName) == 3);
    couplingInterface.setMeshVertices(meshName, 2, pos, vids);
  }
  if (context.isNamed("SolverTwo")) {
    auto meshName = "Test-Square";
    BOOST_REQUIRE(couplingInterface.hasMesh(meshName));
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName) == 3);
    couplingInterface.setMeshVertices(meshName, 2, pos, vids);
  }

  couplingInterface.initialize();
  double dt = couplingInterface.getMaxTimeStepSize();
  while (couplingInterface.isCouplingOngoing()) {
    time += dt;
    couplingInterface.advance(dt);
    dt = couplingInterface.getMaxTimeStepSize();
    timesteps++;
  }
  couplingInterface.finalize();

  BOOST_TEST(time == 10.0);
  BOOST_TEST(timesteps == 10);
}

#endif
