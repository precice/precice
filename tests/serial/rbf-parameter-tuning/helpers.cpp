#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include <boost/proto/proto_fwd.hpp>
#include <iomanip>
#include "mesh/Utils.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/precice.hpp"

std::tuple<std::array<int, 16>, std::array<Eigen::Vector3d, 16>> generateMeshOne(Participant &interface, const std::string &meshOneID)
{
  constexpr double z = 0.3;

  std::array mesh = {
      Eigen::Vector3d{0.0, 0.0, z},
      Eigen::Vector3d{0.5, 0.0, z},
      Eigen::Vector3d{1.0, 0.0, z},
      Eigen::Vector3d{1.5, 0.0, z},
      Eigen::Vector3d{0.0, 0.5, z},
      Eigen::Vector3d{0.5, 0.5, z},
      Eigen::Vector3d{1.0, 0.5, z},
      Eigen::Vector3d{1.5, 0.5, z},
      Eigen::Vector3d{0.0, 1.0, z},
      Eigen::Vector3d{0.5, 1.0, z},
      Eigen::Vector3d{1.0, 1.0, z},
      Eigen::Vector3d{1.5, 1.0, z},
      Eigen::Vector3d{0.0, 1.5, z},
      Eigen::Vector3d{0.5, 1.5, z},
      Eigen::Vector3d{1.0, 1.5, z},
      Eigen::Vector3d{1.5, 1.5, z},
  };

  std::array<int, 16> ids;
  for (size_t i = 0; i < 16; i++) {
    ids[i] = interface.setMeshVertex(meshOneID, mesh[i]);
  }
  return {ids, mesh};
}

std::tuple<std::array<int, 4>, std::array<Eigen::Vector3d, 4>> generateMeshTwo(Participant &interface, const std::string &meshTwoID)
{
  constexpr double z = 0.3;

  std::array mesh = {
      Eigen::Vector3d{0.2, 0.5, z},
      Eigen::Vector3d{1.1, 0.6, z},
      Eigen::Vector3d{0.6, 1.1, z},
      Eigen::Vector3d{1.4, 1.1, z},
  };

  std::array<int, 4> ids;
  for (int i = 0; i < 4; i++) {
    ids[i] = interface.setMeshVertex(meshTwoID, mesh[i]);
  }
  return {ids, mesh};
}

template <size_t N>
std::array<double, N * 3> evaluateFunction(const std::array<Eigen::Vector3d, N> &mesh, double t)
{
  std::array<double, N * 3> values;

  for (size_t i = 0; i < N * 3; i += 3) {
    values[i + 0] = 0.1 * mesh[i / 3].array().sum() + 0;
    values[i + 1] = 0.1 * mesh[i / 3].array().sum() + 1;
    values[i + 2] = 0.1 * mesh[i / 3].array().sum() + 2;
  }

  return values;
}

void testRBFTuning(const std::string configFile, const TestContext &context)
{
  double tA = 0;

  constexpr size_t dim = 3;
  constexpr size_t n   = 4;

  auto meshAID = "MeshOne";
  auto dataAID = "DataOne";
  auto meshBID = "MeshTwo";

  if (context.isNamed("SolverOne")) {
    Participant interfaceA("SolverOne", configFile, 0, 1);
    auto [idsA, meshA] = generateMeshOne(interfaceA, meshAID);
    interfaceA.initialize();

    while (interfaceA.isCouplingOngoing()) {
      double dt = interfaceA.getMaxTimeStepSize();
      interfaceA.writeData(meshAID, dataAID, idsA, evaluateFunction(meshA, tA));
      tA += dt;
      interfaceA.advance(dt);
    }
    interfaceA.finalize();

  } else {
    Participant interfaceB("SolverTwo", configFile, 0, 1);
    auto [idsB, meshB] = generateMeshTwo(interfaceB, meshBID);
    interfaceB.initialize();

    std::array<double, n * dim> expectedValues = evaluateFunction(meshB, tA);

    while (interfaceB.isCouplingOngoing()) {
      double dt = interfaceB.getMaxTimeStepSize();

      std::array<double, n * dim> values;
      for (size_t i = 0; i < values.size(); i++) {
        values[i] = 0;
      }
      interfaceB.readData(meshBID, dataAID, idsB, dt, values);
      interfaceB.advance(dt);

      for (size_t i = 0; i < values.size(); i++) {
        BOOST_TEST(values[i] == expectedValues[i], boost::test_tools::tolerance(1e-5));
      }
    }
    interfaceB.finalize();
  }
}

#endif
