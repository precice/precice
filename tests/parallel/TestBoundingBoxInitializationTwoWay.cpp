#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(TestBoundingBoxInitializationTwoWay)
{
  PRECICE_TEST("Fluid"_on(2_ranks), "Structure"_on(2_ranks));

  std::vector<Eigen::Vector3d> positions;
  std::vector<Eigen::Vector3d> data;
  std::vector<Eigen::Vector3d> expectedData;

  Eigen::Vector3d position;
  Eigen::Vector3d datum;

  int i1 = -1, i2 = -1; //indices for data and positions

  if (context.isNamed("Fluid")) {
    if (context.isPrimary()) {
      i1 = 0;
      i2 = 2;
    } else {
      i1 = 2;
      i2 = 4;
    }
    for (int i = 0; i < 4; i++) {
      position[0] = i * 1.0;
      position[1] = 0;
      position[2] = 0;
      positions.push_back(position);
      datum[0] = i * 1.0;
      datum[1] = i * 2.0;
      datum[2] = i * 3.0;
      data.push_back(datum);
      datum[0] = -i * 1.0;
      datum[1] = -i * 2.0;
      datum[2] = -i * 3.0;
      expectedData.push_back(datum);
    }
  }

  if (context.isNamed("Structure")) {
    if (context.isPrimary()) {
      i1 = 2;
      i2 = 4;
    } else {
      i1 = 0;
      i2 = 2;
    }

    for (int i = 0; i < 4; i++) {
      position[0] = i * 1.0;
      position[1] = 0.0;
      position[2] = 0.0;
      positions.push_back(position);
      datum[0] = -1.0;
      datum[1] = -1.0;
      datum[2] = -1.0;
      data.push_back(datum);
      datum[0] = i * 1.0;
      datum[1] = i * 2.0;
      datum[2] = i * 3.0;
      expectedData.push_back(datum);
    }
  }

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  auto meshName     = context.name + "Mesh";
  auto forcesID     = "Forces";
  auto velocitiesID = "Velocities";

  std::vector<int> vertexIDs;
  for (int i = i1; i < i2; i++) {
    precice::VertexID vertexID = interface.setMeshVertex(meshName, positions[i]);
    vertexIDs.push_back(vertexID);
  }

  interface.initialize();

  if (context.isNamed("Fluid")) {
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      interface.writeData(meshName, forcesID, {&vertexIDs[i], 1}, {data[i + i1].data(), 3});
    }
  }

  double preciceDt = interface.getMaxTimeStepSize();

  if (context.isNamed("Structure")) {
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      interface.readData(meshName, forcesID, {&vertexIDs[i], 1}, preciceDt, {data[i + i1].data(), 3});
      for (size_t d = 0; d < 3; d++) {
        BOOST_TEST(expectedData[i + i1][d] == data[i + i1][d]);
      }
    }

    for (size_t j = 0; j < 4; j++) {
      data[j] = -data[j].array();
    }

    for (size_t i = 0; i < vertexIDs.size(); i++) {
      interface.writeData(meshName, velocitiesID, {&vertexIDs[i], 1}, {data[i + i1].data(), 3});
    }
  }

  interface.advance(1.0);
  preciceDt = interface.getMaxTimeStepSize();

  if (context.isNamed("Fluid")) {
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      interface.readData(meshName, velocitiesID, {&vertexIDs[i], 1}, preciceDt, {data[i + i1].data(), 3});
      for (size_t d = 0; d < 3; d++) {
        BOOST_TEST(expectedData[i + i1][d] == data[i + i1][d]);
      }
    }
  }

  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI
