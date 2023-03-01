#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "precice/SolverInterface.hpp"
#include "testing/Testing.hpp"

/// Tests various distributed communication schemes.
void runTestDistributedCommunication(std::string const &config, TestContext const &context)
{
  std::string meshName;
  int         i1 = -1, i2 = -1; //indices for data and positions

  std::vector<Eigen::VectorXd> positions;
  std::vector<Eigen::VectorXd> data;
  std::vector<Eigen::VectorXd> expectedData;

  Eigen::Vector3d position;
  Eigen::Vector3d datum;

  for (int i = 0; i < 4; i++) {
    position[0] = i * 1.0;
    position[1] = 0.0;
    position[2] = 0.0;
    positions.push_back(position);
    datum[0] = i * 1.0;
    datum[1] = i * 1.0;
    datum[2] = 0.0;
    data.push_back(datum);
    datum[0] = i * 2.0 + 1.0;
    datum[1] = i * 2.0 + 1.0;
    datum[2] = 1.0;
    expectedData.push_back(datum);
  }

  if (context.isNamed("Fluid")) {
    meshName = "FluidMesh";
    if (context.isPrimary()) {
      i1 = 0;
      i2 = 2;
    } else {
      i1 = 2;
      i2 = 4;
    }
  } else {
    meshName = "StructureMesh";
    if (context.isPrimary()) {
      i1 = 0;
      i2 = 1;
    } else {
      i1 = 1;
      i2 = 4;
    }
  }

  precice::SolverInterface precice(context.name, config, context.rank, context.size);
  auto                     meshID   = meshName;
  auto                     forcesID = "Forces";     //  meshID
  auto                     velocID  = "Velocities"; //  meshID

  std::vector<int> vertexIDs;
  for (int i = i1; i < i2; i++) {
    VertexID vertexID = precice.setMeshVertex(meshID, positions[i].data());
    vertexIDs.push_back(vertexID);
  }

  precice.initialize();

  if (context.isNamed("Fluid")) { //Fluid
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.writeVectorData(forcesID, vertexIDs[i], data[i + i1].data());
    }
  } else {
    BOOST_TEST(context.isNamed("Structure"));
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.readVectorData(forcesID, vertexIDs[i], data[i].data());
      data[i] = (data[i] * 2).array() + 1.0;
      precice.writeVectorData(velocID, vertexIDs[i], data[i].data());
    }
  }

  precice.advance(1.0);

  if (context.isNamed("Fluid")) { //Fluid
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      precice.readVectorData(velocID, vertexIDs[i], data[i + i1].data());
      for (size_t d = 0; d < 3; d++) {
        BOOST_TEST(expectedData[i + i1][d] == data[i + i1][d]);
      }
    }
  }

  precice.finalize();
}

#endif
