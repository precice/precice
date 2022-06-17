#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
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

  precice::SolverInterface couplingInterface(context.name, context.config(), 0, 1);
  BOOST_TEST(couplingInterface.getDimensions() == 3);
  if (context.isNamed("SolverOne")) {
    //was necessary to replace pre-defined geometries
    precice::MeshID meshID = couplingInterface.getMeshID("MeshOne");
    couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());

    double dt = couplingInterface.initialize();
    couplingInterface.initializeData();
    while (couplingInterface.isCouplingOngoing()) {
      time += dt;
      dt = couplingInterface.advance(dt);
      timesteps++;
    }
    couplingInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::MeshID meshID = couplingInterface.getMeshID("SolverGeometry");
    int             i0     = couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
    int             i1     = couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(1.0, 0.0, 0.0).data());
    int             i2     = couplingInterface.setMeshVertex(meshID, Eigen::Vector3d(0.0, 1.0, 0.0).data());
    int             e0     = couplingInterface.setMeshEdge(meshID, i0, i1);
    int             e1     = couplingInterface.setMeshEdge(meshID, i1, i2);
    int             e2     = couplingInterface.setMeshEdge(meshID, i2, i0);
    couplingInterface.setMeshTriangle(meshID, e0, e1, e2);
    double dt = couplingInterface.initialize();
    couplingInterface.initializeData();

    int size = couplingInterface.getMeshVertexSize(meshID);
    BOOST_TEST(size == 3);

    while (couplingInterface.isCouplingOngoing()) {
      time += dt;
      dt = couplingInterface.advance(dt);
      timesteps++;
    }
    couplingInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
