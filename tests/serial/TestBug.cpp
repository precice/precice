#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

/**
 * @brief Buggy simulation setup of FSI coupling between Flite and Calculix.
 *
 * Bug: after first call of advance by Flite the mapped forces are value NaN.
 *
 * Some information about the coupling:
 * - explicit coupling scheme
 * - Flite (incompressible Navier-Stokes) starts simulation
 * - Mapping is done on Flite side with RBF
 *
 * @todo rename this test and config
 */
BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(TestBug)
{
  PRECICE_TEST("Flite"_on(1_rank), "Calculix"_on(1_rank));
  using Eigen::Vector3d;

  int                   slices = 5;
  std::vector<Vector3d> coords;
  for (int i = 0; i < slices; i++) {
    double z = (double) i * 1.0;
    coords.push_back(Vector3d(1.0, 0.0, z));
    coords.push_back(Vector3d(0.0, 1.0, z));
    coords.push_back(Vector3d(-1.0, 0.0, z));
    coords.push_back(Vector3d(0.0, -1.0, z));
  }

  if (context.isNamed("Flite")) {
    precice::SolverInterface interface(context.name, context.config(), 0, 1);

    precice::MeshID meshID             = interface.getMeshID("FliteNodes");
    int             forcesID           = interface.getDataID("Forces", meshID);
    int             displacementsID    = interface.getDataID("Displacements", meshID);
    int             oldDisplacementsID = interface.getDataID("OldDisplacements", meshID);
    BOOST_TEST(interface.getDimensions() == 3);
    for (Vector3d &coord : coords) {
      interface.setMeshVertex(meshID, coord.data());
    }
    double maxDt = interface.initialize();
    double dt    = 1.0e-5 / 15.0; // Flite took 15 subcycling steps
    while (interface.isCouplingOngoing()) {
      dt = dt < maxDt ? dt : maxDt;
      for (int i = 0; i < (int) coords.size(); i++) {
        double force[3] = {1.0, 2.0, 3.0};
        interface.writeVectorData(forcesID, i, force);
      }
      maxDt = interface.advance(dt);
      interface.mapReadDataTo(meshID);
      for (int i = 0; i < (int) coords.size(); i++) {
        double displacement[3];
        double oldDisplacement[3];
        interface.readVectorData(displacementsID, i, displacement);
        interface.readVectorData(oldDisplacementsID, i, oldDisplacement);
      }
    }
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("Calculix"));
    precice::SolverInterface interface(context.name, context.config(), 0, 1);

    precice::MeshID meshID = interface.getMeshID("CalculixNodes");
    for (Vector3d &coord : coords) {
      interface.setMeshVertex(meshID, coord.data());
    }
    for (int i = 0; i < slices - 1; i++) {
      // Build cylinder/channel geometry
      interface.setMeshTriangleWithEdges(meshID, i * 4, (i * 4) + 1, (i + 1) * 4);
      interface.setMeshTriangleWithEdges(meshID, (i + 1) * 4, (i * 4) + 1, ((i + 1) * 4) + 1);
      interface.setMeshTriangleWithEdges(meshID, i * 4 + 1, (i * 4) + 2, (i + 1) * 4 + 1);
      interface.setMeshTriangleWithEdges(meshID, (i + 1) * 4 + 1, (i * 4) + 2, ((i + 1) * 4) + 2);
      interface.setMeshTriangleWithEdges(meshID, i * 4 + 2, (i * 4) + 3, (i + 1) * 4 + 2);
      interface.setMeshTriangleWithEdges(meshID, (i + 1) * 4 + 2, (i * 4) + 3, ((i + 1) * 4) + 3);
      interface.setMeshTriangleWithEdges(meshID, i * 4 + 3, (i * 4), (i + 1) * 4 + 3);
      interface.setMeshTriangleWithEdges(meshID, (i + 1) * 4 + 3, i * 4, (i + 1) * 4);
    }
    double dt = interface.initialize();
    while (interface.isCouplingOngoing()) {
      interface.advance(dt);
    }
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
