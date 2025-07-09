#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "testing/Testing.hpp"

#include "mesh/Utils.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/precice.hpp"

std::vector<int> generateMeshOne(precice::Participant &interface, const std::string &meshOneID)
{
  const double     z = 0.3;
  std::vector<int> ids;
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{0.0, 0.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{1.0, 0.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{1.0, 1.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{0.0, 1.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{2.0, 0.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{3.0, 0.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{3.0, 1.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{2.0, 1.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{4.0, 0.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{5.0, 0.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{5.0, 1.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshOneID, Eigen::Vector3d{4.0, 1.0, z}));
  return ids;
}

std::vector<int> generateMeshTwo(precice::Participant &interface, const std::string &meshTwoID)
{
  const double     z = 0.3;
  std::vector<int> ids;
  ids.emplace_back(interface.setMeshVertex(meshTwoID, Eigen::Vector3d{0.0, 0.0, z}));
  ids.emplace_back(interface.setMeshVertex(meshTwoID, Eigen::Vector3d{0.5, 0.5, z}));
  ids.emplace_back(interface.setMeshVertex(meshTwoID, Eigen::Vector3d{3.5, 0.5, z}));
  return ids;
}


std::vector<double> evaluateFunction(double t) {
  std::vector<double> values;
  double a = 0.5 + 0.5 * std::cos(2*t);

  for (unsigned int i = 0; i < 12; ++i) {
    values.emplace_back(a * std::pow(i + 1, 2));
    values.emplace_back(a * i + 1);
    values.emplace_back(a);
  }
  return values;
}


void testRBFTuning(const std::string configFile, const TestContext &context)
{
  std::array<double, 5 * 9> expectedValues = {
    1.0000000000000002,  1,                  1,                   2.187913147209069,   0.7361819229699296,  0.2957923081147373,   22.859664317930537, 2.507197240750358,  0.29496764858906593,
    0.2919265817264288,  1,                  0.29192658172642877, 0.6387100061790564,  0.424353743007232,   0.08634963740890586,  6.6733436637470795, 0.9407762714227634, 0.08610889737248849,
    0.17317818956819403, 0.9999999999999999, 0.17317818956819406, 0.37889883776611616, 0.3720581843199937,  0.051224776407507625, 3.958795280715856,  0.6780775642687428, 0.0510819633638417,
    0.9800851433251829,  1,                  0.9800851433251829,  2.1443411704654523,  0.7274116269090105,  0.2899016466931189,   22.404417379404524, 2.4631410054906766, 0.2890934101437068,
    0.42724998309569323, 1,                  0.42724998309569323, 0.9347858551599195,  0.48394876361713696, 0.12637725864185756,  9.766791193409045,  1.24014270444377,   0.1260249228734548
  };

  double tA = 0;

  auto meshAID = "MeshOne";
  auto dataAID = "DataOne";
  auto meshBID = "MeshTwo";

  if (context.isNamed("SolverOne")) {
    precice::Participant interfaceA("SolverOne", configFile, 0, 1);
    std::vector<int> idsA = generateMeshOne(interfaceA, meshAID);
    interfaceA.initialize();

    while (interfaceA.isCouplingOngoing()) {
      double dt = interfaceA.getMaxTimeStepSize();
      interfaceA.writeData(meshAID, dataAID, idsA, evaluateFunction(tA));
      tA += dt;
      interfaceA.advance(dt);
    }
    interfaceA.finalize();
  } else {
    precice::Participant interfaceB("SolverTwo", configFile, 0, 1);
    std::vector<int> idsB = generateMeshTwo(interfaceB, meshBID);
    interfaceB.initialize();

    int it = 0;

    while (interfaceB.isCouplingOngoing()) {

      double dt = interfaceB.getMaxTimeStepSize();
      std::array<double, 9> values;
      interfaceB.readData(meshBID, dataAID, idsB, dt, values);
      interfaceB.advance(dt);

      for (size_t i = 0; i < values.size(); i++) {
        BOOST_TEST(values[i] == expectedValues[i + it * 9], boost::test_tools::tolerance(1e-7));
      }
      it++;
    }
    interfaceB.finalize();
  }
}

#endif