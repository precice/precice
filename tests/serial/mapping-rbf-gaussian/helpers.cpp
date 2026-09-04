#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include <boost/proto/proto_fwd.hpp>
#include <iomanip>
#include "mesh/Utils.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/precice.hpp"

std::tuple<std::array<int, 12>, std::array<Eigen::Vector3d, 12>> generateMeshOneGauss(Participant &interface, const std::string &meshOneID)
{
  constexpr int    N = 12;
  constexpr double z = 0.3;

  std::array mesh = {
      Eigen::Vector3d{0.0, 0.0, z},
      Eigen::Vector3d{1.0, 0.0, z},
      Eigen::Vector3d{1.0, 1.0, z},
      Eigen::Vector3d{0.0, 1.0, z},
      Eigen::Vector3d{2.0, 0.0, z},
      Eigen::Vector3d{3.0, 0.0, z},
      Eigen::Vector3d{3.0, 1.0, z},
      Eigen::Vector3d{2.0, 1.0, z},
      Eigen::Vector3d{4.0, 0.0, z},
      Eigen::Vector3d{5.0, 0.0, z},
      Eigen::Vector3d{5.0, 1.0, z},
      Eigen::Vector3d{4.0, 1.0, z},
  };

  std::array<int, N> ids;
  for (size_t i = 0; i < N; i++) {
    ids[i] = interface.setMeshVertex(meshOneID, mesh[i]);
  }
  return {ids, mesh};
}

std::tuple<std::array<int, 3>, std::array<Eigen::Vector3d, 3>> generateMeshTwoGauss(Participant &interface, const std::string &meshTwoID)
{
  constexpr int    N = 3;
  constexpr double z = 0.3;

  std::array mesh = {
      Eigen::Vector3d{0.0, 0.0, z}, // Maps to vertex A
      Eigen::Vector3d{0.5, 0.5, z}, // Maps on the left side of the domain
      Eigen::Vector3d{3.5, 0.5, z}, // Maps more in the middle of the domain
  };

  std::array<int, N> ids;
  for (int i = 0; i < N; i++) {
    ids[i] = interface.setMeshVertex(meshTwoID, mesh[i]);
  }
  return {ids, mesh};
}

template <size_t N>
std::array<double, N * 3> evaluateFunctionGauss(const std::array<Eigen::Vector3d, N> &mesh)
{
  std::array<double, N * 3> values;

  for (size_t i = 0; i < N * 3; i += 3) {
    values[i + 0] = 0.1 * mesh[i / 3].array().sum() + 0;
    values[i + 1] = 0.1 * mesh[i / 3].array().sum() + 1;
    values[i + 2] = 0.1 * mesh[i / 3].array().sum() + 2;
  }
  return values;
}

void testRBFMappingVectorial(const std::string configFile, const TestContext &context, bool mappingIsConservative)
{
  constexpr size_t dim    = 3;
  constexpr size_t numOut = 3;
  constexpr size_t numIn  = 12;

  auto meshOneID = "MeshOne";
  auto dataOneID = "DataOne";
  auto meshTwoID = "MeshTwo";

  if (context.isNamed("SolverOne")) {
    Participant interfaceA("SolverOne", configFile, 0, 1);

    auto [idsOne, meshOne] = generateMeshOneGauss(interfaceA, meshOneID);
    interfaceA.initialize();

    while (interfaceA.isCouplingOngoing()) {
      double maxDt = interfaceA.getMaxTimeStepSize();

      std::array<double, numIn * dim> inputValues = evaluateFunctionGauss(meshOne);
      interfaceA.writeData(meshOneID, dataOneID, idsOne, inputValues);
      interfaceA.advance(maxDt);
    }
    interfaceA.finalize();

  } else {
    Participant interfaceB("SolverTwo", configFile, 0, 1);

    auto [idsTwo, meshTwo] = generateMeshTwoGauss(interfaceB, meshTwoID);
    interfaceB.initialize();

    std::array<double, numOut * dim> expectedValues;
    if (mappingIsConservative) {
      expectedValues = {-0.6, -0.6, -0.6, 0.853333, 4.85333, 8.85333, 3.70667, 11.7067, 19.7067};
    } else {
      expectedValues = evaluateFunctionGauss(meshTwo);
    }

    double maxDt = interfaceB.getMaxTimeStepSize();

    Eigen::Vector<double, numOut * dim> values = Eigen::Vector<double, numOut * dim>::Zero();
    interfaceB.readData(meshTwoID, dataOneID, idsTwo, maxDt, values);
    interfaceB.advance(maxDt);

    for (Eigen::Index i = 0; i < values.size(); i++) {
      BOOST_TEST(values[i] == expectedValues[i], boost::test_tools::tolerance(1e-5));
    }
    interfaceB.finalize();
  }
}

void testRBFMapping(const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  std::vector<double> values;
  const unsigned int  nCoords = 12;
  for (unsigned int i = 0; i < nCoords; ++i)
    values.emplace_back(std::pow(i + 1, 2));

  double expectedValTwoA = 1.0000000014191541;
  double expectedValTwoB = 7.30892688709867;
  double expectedValTwoC = 77.5938805368033;

  if (context.isNamed("SolverOne")) {
    precice::Participant interface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshOneID = "MeshOne";

    // Setup mesh one.
    auto [ids, _] = generateMeshOneGauss(interface, meshOneID);

    // Initialize, thus sending the mesh.
    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant should have to advance once!");
    // Write the data to be send.
    auto dataAID = "DataOne";
    BOOST_TEST(!interface.requiresGradientDataFor(meshOneID, dataAID));

    interface.writeData(meshOneID, dataAID, ids, values);

    // Advance, thus send the data to the receiving partner.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant should have to advance once!");
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::Participant interface("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshTwoID = "MeshTwo";

    // Setup receiving mesh.
    auto [ids, _] = generateMeshTwoGauss(interface, meshTwoID);

    // Initialize, thus receive the data and map.
    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    auto dataAID = "DataOne";
    BOOST_TEST(!interface.requiresGradientDataFor(meshTwoID, dataAID));

    double values[3];
    interface.readData(meshTwoID, dataAID, ids, maxDt, values);

    // Due to Eigen 3.3.7 (Ubunu 2004) giving slightly different results
    BOOST_TEST(values[0] == expectedValTwoA, boost::test_tools::tolerance(1e-8));
    BOOST_TEST(values[1] == expectedValTwoB, boost::test_tools::tolerance(3e-2));
    BOOST_TEST(values[2] == expectedValTwoC, boost::test_tools::tolerance(2e-3));

    // Verify that there is only one time step necessary.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant should have to advance once!");
    interface.finalize();
  }
}

#endif
