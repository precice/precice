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

    int idA = interface.setMeshVertex(meshName, coordA);
    int idB = interface.setMeshVertex(meshName, coordB);
    int idC = interface.setMeshVertex(meshName, coordC);
    int idD = interface.setMeshVertex(meshName, coordD);

    // Initialize, the mesh
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Read the summed data from the mesh.
    auto   dataAID = "Target";
    double valueA, valueB, valueC, valueD;

    while (interface.isCouplingOngoing()) {

      interface.readData(meshName, dataAID, {&idA, 1}, dt, {&valueA, 1});
      interface.readData(meshName, dataAID, {&idB, 1}, dt, {&valueB, 1});
      interface.readData(meshName, dataAID, {&idC, 1}, dt, {&valueC, 1});
      interface.readData(meshName, dataAID, {&idD, 1}, dt, {&valueD, 1});

      BOOST_TEST(valueA == expectedValueA);
      BOOST_TEST(valueB == expectedValueB);
      BOOST_TEST(valueC == expectedValueC);
      BOOST_TEST(valueD == expectedValueD);

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
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

    int idA = interface.setMeshVertex(meshName, coordA);
    int idB = interface.setMeshVertex(meshName, coordB);
    int idC = interface.setMeshVertex(meshName, coordC);
    int idD = interface.setMeshVertex(meshName, coordD);

    // Initialize, the mesh
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    auto   dataAID = "SourceOne";
    double valueA  = 1.0;
    double valueB  = 3.0;
    double valueC  = 5.0;
    double valueD  = 7.0;

    while (interface.isCouplingOngoing()) {

      interface.writeData(meshName, dataAID, {&idA, 1}, {&valueA, 1});
      interface.writeData(meshName, dataAID, {&idB, 1}, {&valueB, 1});
      interface.writeData(meshName, dataAID, {&idC, 1}, {&valueC, 1});
      interface.writeData(meshName, dataAID, {&idD, 1}, {&valueD, 1});

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
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

    int idA = interface.setMeshVertex(meshName, coordA);
    int idB = interface.setMeshVertex(meshName, coordB);
    int idC = interface.setMeshVertex(meshName, coordC);
    int idD = interface.setMeshVertex(meshName, coordD);

    // Initialize, the mesh
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    auto   dataAID = "SourceTwo";
    double valueA  = 2.0;
    double valueB  = 4.0;
    double valueC  = 6.0;
    double valueD  = 8.0;

    while (interface.isCouplingOngoing()) {

      interface.writeData(meshName, dataAID, {&idA, 1}, {&valueA, 1});
      interface.writeData(meshName, dataAID, {&idB, 1}, {&valueB, 1});
      interface.writeData(meshName, dataAID, {&idC, 1}, {&valueC, 1});
      interface.writeData(meshName, dataAID, {&idD, 1}, {&valueD, 1});

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
    }

    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
