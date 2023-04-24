#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(SummationActionTwoSources)
{
  PRECICE_TEST("SolverTarget"_on(1_rank), "SolverSourceOne"_on(1_rank), "SolverSourceTwo"_on(1_rank));

  using Eigen::Vector3d;

  if (context.isNamed("SolverTarget")) {
    // Expected values in the target solver
    double expectedValueA = 3.0;
    double expectedValueB = 7.0;
    double expectedValueC = 11.0;
    double expectedValueD = 15.0;

    // Target solver
    precice::SolverInterface interface(context.name, context.config(), 0, 1);

    // Set mesh
    Vector3d coordA{0.0, 0.0, 0.3};
    Vector3d coordB{1.0, 0.0, 0.3};
    Vector3d coordC{1.0, 1.0, 0.3};
    Vector3d coordD{0.0, 1.0, 0.3};

    auto meshName = "MeshTarget";

    int idA = interface.setMeshVertex(meshName, coordA.data());
    int idB = interface.setMeshVertex(meshName, coordB.data());
    int idC = interface.setMeshVertex(meshName, coordC.data());
    int idD = interface.setMeshVertex(meshName, coordD.data());

    // Initialize, the mesh
    double dt = interface.initialize();

    // Read the summed data from the mesh.
    auto   dataAID = "Target";
    double valueA, valueB, valueC, valueD;

    while (interface.isCouplingOngoing()) {

      interface.readScalarData(meshName, dataAID, idA, dt, valueA);
      interface.readScalarData(meshName, dataAID, idB, dt, valueB);
      interface.readScalarData(meshName, dataAID, idC, dt, valueC);
      interface.readScalarData(meshName, dataAID, idD, dt, valueD);

      BOOST_TEST(valueA == expectedValueA);
      BOOST_TEST(valueB == expectedValueB);
      BOOST_TEST(valueC == expectedValueC);
      BOOST_TEST(valueD == expectedValueD);

      dt = interface.advance(dt);
    }

    interface.finalize();
  } else if (context.isNamed("SolverSourceOne")) {
    // Source solver one
    precice::SolverInterface interface(context.name, context.config(), 0, 1);

    // Set mesh
    Vector3d coordA{0.0, 0.0, 0.3};
    Vector3d coordB{1.0, 0.0, 0.3};
    Vector3d coordC{1.0, 1.0, 0.3};
    Vector3d coordD{0.0, 1.0, 0.3};

    auto meshName = "MeshOne";

    int idA = interface.setMeshVertex(meshName, coordA.data());
    int idB = interface.setMeshVertex(meshName, coordB.data());
    int idC = interface.setMeshVertex(meshName, coordC.data());
    int idD = interface.setMeshVertex(meshName, coordD.data());

    // Initialize, the mesh
    double dt = interface.initialize();

    auto   dataAID = "SourceOne";
    double valueA  = 1.0;
    double valueB  = 3.0;
    double valueC  = 5.0;
    double valueD  = 7.0;

    while (interface.isCouplingOngoing()) {

      interface.writeScalarData(meshName, dataAID, idA, valueA);
      interface.writeScalarData(meshName, dataAID, idB, valueB);
      interface.writeScalarData(meshName, dataAID, idC, valueC);
      interface.writeScalarData(meshName, dataAID, idD, valueD);

      dt = interface.advance(dt);
    }
    interface.finalize();
  } else {
    BOOST_REQUIRE(context.isNamed("SolverSourceTwo"));
    // Source solver two
    precice::SolverInterface interface(context.name, context.config(), 0, 1);

    // Set mesh
    Vector3d coordA{0.0, 0.0, 0.3};
    Vector3d coordB{1.0, 0.0, 0.3};
    Vector3d coordC{1.0, 1.0, 0.3};
    Vector3d coordD{0.0, 1.0, 0.3};

    auto meshName = "MeshTwo";

    int idA = interface.setMeshVertex(meshName, coordA.data());
    int idB = interface.setMeshVertex(meshName, coordB.data());
    int idC = interface.setMeshVertex(meshName, coordC.data());
    int idD = interface.setMeshVertex(meshName, coordD.data());

    // Initialize, the mesh
    double dt = interface.initialize();

    auto   dataAID = "SourceTwo";
    double valueA  = 2.0;
    double valueB  = 4.0;
    double valueC  = 6.0;
    double valueD  = 8.0;

    while (interface.isCouplingOngoing()) {

      interface.writeScalarData(meshName, dataAID, idA, valueA);
      interface.writeScalarData(meshName, dataAID, idB, valueB);
      interface.writeScalarData(meshName, dataAID, idC, valueC);
      interface.writeScalarData(meshName, dataAID, idD, valueD);

      dt = interface.advance(dt);
    }

    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
