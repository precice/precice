#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
/**
  * @brief Runs a coupled simulation where one solver supplies a geometry.
  *
  * SolverOne only reads the displacements of the geometry and checks whether
  * they are equals to the coordinates of SolverTwo. SolverTwo creates and
  * displaces the coordinates.
  *
  * @todo Maybe remove this test.
  */
BOOST_AUTO_TEST_CASE(TestExplicitWithSolverGeometry)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  int    timesteps = 0;
  double time      = 0;

  double v0[]   = {0, 0, 0};
  double v100[] = {1, 0, 0};
  double v010[] = {0, 1, 0};

  precice::Participant couplingInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    //was necessary to replace pre-defined geometries
    auto meshName = "MeshOne";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName) == 3);
    couplingInterface.setMeshVertex(meshName, v0);
    couplingInterface.setMeshVertex(meshName, v100);

    couplingInterface.initialize();
    double dt = couplingInterface.getMaxTimeStepSize();
    while (couplingInterface.isCouplingOngoing()) {
      time += dt;
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
      timesteps++;
    }
    couplingInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto meshName = "SolverGeometry";
    int  i0       = couplingInterface.setMeshVertex(meshName, v0);
    int  i1       = couplingInterface.setMeshVertex(meshName, v100);
    int  i2       = couplingInterface.setMeshVertex(meshName, v010);
    couplingInterface.setMeshTriangle(meshName, i0, i1, i2);
    couplingInterface.initialize();
    double dt = couplingInterface.getMaxTimeStepSize();

    int size = couplingInterface.getMeshVertexSize(meshName);
    BOOST_TEST(size == 3);

    while (couplingInterface.isCouplingOngoing()) {
      time += dt;
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
      timesteps++;
    }
    couplingInterface.finalize();
    BOOST_TEST(time == 0.05);
    BOOST_TEST(timesteps == 5);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
