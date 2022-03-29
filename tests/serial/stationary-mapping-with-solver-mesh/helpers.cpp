#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "precice/SolverInterface.hpp"
#include "testing/Testing.hpp"

void runTestStationaryMappingWithSolverMesh(std::string const &config, int dim, TestContext const &context)
{
  std::string meshForcesA = "MeshForcesA";
  std::string meshDisplA  = "MeshDisplacementsA";
  std::string meshForcesB = "MeshForcesB";
  std::string meshDisplB  = "MeshDisplacementsB";
  std::string dataForces  = "Forces";
  std::string dataDispl   = "Displacements";
  using testing::equals;

  precice::SolverInterface interface(context.name, config, 0, 1);
  BOOST_TEST(interface.getDimensions() == dim);

  std::vector<Eigen::VectorXd> positions;
  Eigen::VectorXd              position(dim);
  if (dim == 2) {
    position << 0.0, 0.0;
    positions.push_back(position);
    position << 1.0, 0.0;
    positions.push_back(position);
    position << 1.0, 1.0;
    positions.push_back(position);
    position << 0.0, 1.0;
    positions.push_back(position);
  } else {
    position << 0.0, 0.0, 0.0;
    positions.push_back(position);
    position << 1.0, 0.0, 0.0;
    positions.push_back(position);
    position << 1.0, 1.0, 0.0;
    positions.push_back(position);
    position << 0.0, 1.0, 1.0;
    positions.push_back(position);
    position << 0.0, 0.0, 1.0;
    positions.push_back(position);
  }
  size_t size = positions.size();

  if (context.isNamed("SolverA")) {
    int meshForcesID = interface.getMeshID(meshForcesA);
    int meshDisplID  = interface.getMeshID(meshDisplA);
    int dataForcesID = interface.getDataID(dataForces, meshForcesID);
    int dataDisplID  = interface.getDataID(dataDispl, meshDisplID);

    // Set solver mesh positions for reading and writing data with mappings
    for (size_t i = 0; i < size; i++) {
      position = positions.at(i).array() + 0.1;
      interface.setMeshVertex(meshForcesID, position.data());
      position = positions.at(i).array() + 0.6;
      interface.setMeshVertex(meshDisplID, position.data());
    }
    double maxDt = interface.initialize();

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(not interface.isReadDataAvailable());
    Eigen::VectorXd force = Eigen::VectorXd::Constant(dim, 1);
    Eigen::VectorXd displ = Eigen::VectorXd::Constant(dim, 0);
    for (size_t i = 0; i < size; i++) {
      interface.writeVectorData(dataForcesID, i, force.data());
    }
    interface.mapWriteDataFrom(meshForcesID);
    maxDt = interface.advance(maxDt);
    interface.mapReadDataTo(meshDisplID);

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(interface.isReadDataAvailable());
    force.array() += 1.0;
    for (size_t i = 0; i < size; i++) {
      interface.readVectorData(dataDisplID, i, displ.data());
      BOOST_TEST(displ(0) == positions.at(i)(0) + 0.1);
      interface.writeVectorData(dataForcesID, i, force.data());
    }
    interface.mapWriteDataFrom(meshForcesID);
    maxDt = interface.advance(maxDt);
    interface.mapReadDataTo(meshDisplID);

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(interface.isReadDataAvailable());
    for (size_t i = 0; i < size; i++) {
      interface.readVectorData(dataDisplID, i, displ.data());
      BOOST_TEST(displ(0) == 2.0 * (positions.at(i)(0) + 0.1));
    }
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverB"));
    int meshForcesID = interface.getMeshID(meshForcesB);
    int meshDisplID  = interface.getMeshID(meshDisplB);
    int dataForcesID = interface.getDataID(dataForces, meshForcesID);
    int dataDisplID  = interface.getDataID(dataDispl, meshDisplID);

    // Set solver mesh positions provided to SolverA for data mapping
    for (size_t i = 0; i < size; i++) {
      interface.setMeshVertex(meshForcesID, positions.at(i).data());
      position = positions.at(i).array() + 0.5;
      interface.setMeshVertex(meshDisplID, position.data());
    }
    double maxDt = interface.initialize();

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(interface.isReadDataAvailable());
    Eigen::VectorXd force      = Eigen::VectorXd::Zero(dim);
    Eigen::VectorXd totalForce = Eigen::VectorXd::Zero(dim);
    Eigen::VectorXd displ      = Eigen::VectorXd::Zero(dim);
    for (size_t i = 0; i < size; i++) {
      interface.readVectorData(dataForcesID, i, force.data());
      totalForce += force;
      displ.setConstant(positions.at(i)(0));
      interface.writeVectorData(dataDisplID, i, displ.data());
    }
    Eigen::VectorXd expected = Eigen::VectorXd::Constant(dim, size);
    BOOST_TEST(equals(totalForce, expected));
    maxDt = interface.advance(maxDt);

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(interface.isReadDataAvailable());
    totalForce.setConstant(0);
    for (size_t i = 0; i < positions.size(); i++) {
      interface.readVectorData(dataForcesID, i, force.data());
      totalForce += force;
      displ.setConstant(2.0 * positions.at(i)(0));
      interface.writeVectorData(dataDisplID, i, displ.data());
    }
    expected.setConstant(2.0 * (double) size);
    BOOST_TEST(equals(totalForce, expected));
    maxDt = interface.advance(maxDt);

    BOOST_TEST(interface.isWriteDataRequired(maxDt));
    BOOST_TEST(not interface.isReadDataAvailable()); //second participant has no new data after last advance
    for (size_t i = 0; i < size; i++) {
      interface.readVectorData(dataForcesID, i, force.data());
    }
    interface.finalize();
  }
}

#endif