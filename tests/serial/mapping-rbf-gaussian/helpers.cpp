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
  Vector3d coordOneE{2.0, 0.0, z};
  Vector3d coordOneF{3.0, 0.0, z};
  Vector3d coordOneG{3.0, 1.0, z};
  Vector3d coordOneH{2.0, 1.0, z};
  Vector3d coordOneI{4.0, 0.0, z};
  Vector3d coordOneJ{5.0, 0.0, z};
  Vector3d coordOneK{5.0, 1.0, z};
  Vector3d coordOneL{4.0, 1.0, z};

  std::vector<double> values;
  const unsigned int  nCoords = 12;
  for (uint i = 0; i < nCoords; ++i)
    values.emplace_back(std::pow(i + 1, 2));
  // MeshTwo
  Vector3d coordTwoA{0.0, 0.0, z}; // Maps to vertex A
  Vector3d coordTwoB{0.5, 0.5, z}; // Maps on the left side of the domain
  Vector3d coordTwoC{3.5, 0.5, z}; // Maps more in the middle of the domain

  double expectedValTwoA = 1.0000000014191541;
  double expectedValTwoB = 7.1554050349583402;
  double expectedValTwoC = 77.68217404046861;

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface interface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    const int meshOneID = interface.getMeshID("MeshOne");

    // Setup mesh one.
    std::vector<int> ids;
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneA.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneB.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneC.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneD.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneE.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneF.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneG.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneH.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneI.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneJ.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneK.data()));
    ids.emplace_back(interface.setMeshVertex(meshOneID, coordOneL.data()));

    // Initialize, thus sending the mesh.
    double maxDt = interface.initialize();
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant should have to advance once!");
    // Write the data to be send.
    int dataAID = interface.getDataID("DataOne", meshOneID);
    BOOST_TEST(!interface.isGradientDataRequired(dataAID));

    interface.writeBlockScalarData(dataAID, nCoords, ids.data(), values.data());

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
    int dataAID = interface.getDataID("DataOne", meshTwoID);
    BOOST_TEST(!interface.isGradientDataRequired(dataAID));

    double valueA, valueB, valueC;
    interface.readScalarData(dataAID, idA, valueA);
    interface.readScalarData(dataAID, idB, valueB);
    interface.readScalarData(dataAID, idC, valueC);

    // Due to Eigen 3.3.7 (Ubunu 2004) giving slightly different results
    BOOST_TEST(valueA == expectedValTwoA, boost::test_tools::tolerance(1e-8));
    BOOST_TEST(valueB == expectedValTwoB);
    BOOST_TEST(valueC == expectedValTwoC);

    // Verify that there is only one time step necessary.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant should have to advance once!");
    interface.finalize();
  }
}

#endif
