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
  for (unsigned int i = 0; i < nCoords; ++i)
    values.emplace_back(std::pow(i + 1, 2));
  // MeshTwo
  Vector3d coordTwoA{0.0, 0.0, z}; // Maps to vertex A
  Vector3d coordTwoB{0.5, 0.5, z}; // Maps on the left side of the domain
  Vector3d coordTwoC{3.5, 0.5, z}; // Maps more in the middle of the domain

  double expectedValTwoA = 1.0000000014191541;
  double expectedValTwoB = 7.30892688709867;
  double expectedValTwoC = 77.5938805368033;

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface interface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshOneID = "MeshOne";

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
    auto dataAID = "DataOne"; //  meshOneID
    BOOST_TEST(!interface.requiresGradientDataFor(meshID, dataAID));

    interface.writeBlockScalarData(meshID, dataAID, nCoords, ids.data(), values.data());

    // Advance, thus send the data to the receiving partner.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant should have to advance once!");
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::SolverInterface interface("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshTwoID = "MeshTwo";

    // Setup receiving mesh.
    int idA = interface.setMeshVertex(meshTwoID, coordTwoA.data());
    int idB = interface.setMeshVertex(meshTwoID, coordTwoB.data());
    int idC = interface.setMeshVertex(meshTwoID, coordTwoC.data());

    // Initialize, thus receive the data and map.
    double maxDt = interface.initialize();
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    auto dataAID = "DataOne"; //  meshTwoID
    BOOST_TEST(!interface.requiresGradientDataFor(meshID, dataAID));

    double valueA, valueB, valueC;
    interface.readScalarData(meshID, dataAID, idA, valueA);
    interface.readScalarData(meshID, dataAID, idB, valueB);
    interface.readScalarData(meshID, dataAID, idC, valueC);

    // Due to Eigen 3.3.7 (Ubunu 2004) giving slightly different results
    BOOST_TEST(valueA == expectedValTwoA, boost::test_tools::tolerance(1e-8));
    BOOST_TEST(valueB == expectedValTwoB, boost::test_tools::tolerance(3e-2));
    BOOST_TEST(valueC == expectedValTwoC, boost::test_tools::tolerance(2e-3));

    // Verify that there is only one time step necessary.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant should have to advance once!");
    interface.finalize();
  }
}

#endif
