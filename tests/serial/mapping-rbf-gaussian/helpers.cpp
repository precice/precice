#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "mesh/Utils.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"

void testRBFMapping(const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  const double z = 0.3;

  // MeshOne
  Vector3d coordOneA{0.0, 0.0, z};
  Vector3d coordOneB{1.0, 0.0, z};
  Vector3d coordOneC{1.0, 1.0, z};
  Vector3d coordOneD{0.0, 1.0, z};
  double   valOneA = 1.0;
  double   valOneB = 3.0;
  double   valOneC = 5.0;
  double   valOneD = 7.0;

  // MeshTwo
  Vector3d coordTwoA{0.0, 0.0, z}; // Maps to vertex A
  Vector3d coordTwoB{0.0, 0.5, z}; // Maps to edge AD

  double expectedValTwoA = 1.0;
  double expectedValTwoB = 4.0;

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface interface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    const int meshOneID = interface.getMeshID("MeshOne");

    // Setup mesh one.
    int idA = interface.setMeshVertex(meshOneID, coordOneA.data());
    int idB = interface.setMeshVertex(meshOneID, coordOneB.data());
    int idC = interface.setMeshVertex(meshOneID, coordOneC.data());
    int idD = interface.setMeshVertex(meshOneID, coordOneD.data());

    // Initialize, thus sending the mesh.
    double maxDt = interface.initialize();
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant should have to advance once!");

    // Write the data to be send.
    int dataAID = interface.getDataID("DataOne", meshOneID);
    interface.writeScalarData(dataAID, idA, valOneA);
    interface.writeScalarData(dataAID, idB, valOneB);
    interface.writeScalarData(dataAID, idC, valOneC);
    interface.writeScalarData(dataAID, idD, valOneD);

    // Advance, thus send the data to the receiving partner.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant should have to advance once!");
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::SolverInterface interface("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    int meshTwoID = interface.getMeshID("MeshTwo");

    // Setup receiving mesh.
    int idA = interface.setMeshVertex(meshTwoID, coordTwoA.data());
    int idB = interface.setMeshVertex(meshTwoID, coordTwoB.data());
    int idC = interface.setMeshVertex(meshTwoID, coordTwoC.data());

    // Initialize, thus receive the data and map.
    double maxDt = interface.initialize();
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    int    dataAID = interface.getDataID("DataOne", meshTwoID);
    double valueA, valueB, valueC;
    interface.readScalarData(dataAID, idA, valueA);
    interface.readScalarData(dataAID, idB, valueB);
    interface.readScalarData(dataAID, idC, valueC);

    BOOST_TEST(valueA == expectedValTwoA);
    BOOST_TEST(valueB == expectedValTwoB);

    // Verify that there is only one time step necessary.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant should have to advance once!");
    interface.finalize();
  }
}

#endif