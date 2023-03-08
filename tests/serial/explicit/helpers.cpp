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

  //was necessary to replace pre-defined geometries
  if (context.isNamed("SolverOne") && couplingInterface.hasMesh("MeshOne")) {
    auto meshName = "MeshOne";
    couplingInterface.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    couplingInterface.setMeshVertex(meshName, Eigen::Vector3d(1.0, 0.0, 0.0).data());
  }
  if (context.isNamed("SolverTwo") && couplingInterface.hasMesh("Test-Square")) {
    auto meshName = "Test-Square";
    couplingInterface.setMeshVertex(meshName, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    couplingInterface.setMeshVertex(meshName, Eigen::Vector3d(1.0, 0.0, 0.0).data());
  }

  BOOST_TEST(couplingInterface.getDimensions() == 3);
  double dt = couplingInterface.initialize();
  while (couplingInterface.isCouplingOngoing()) {
    time += dt;
    dt = couplingInterface.advance(dt);
    timesteps++;
  }
  couplingInterface.finalize();

  BOOST_TEST(time == 10.0);
  BOOST_TEST(timesteps == 10);
}

#endif
