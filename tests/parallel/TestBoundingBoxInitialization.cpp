#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
PRECICE_TEST_SETUP("Fluid"_on(2_ranks), "Structure"_on(2_ranks))
BOOST_AUTO_TEST_CASE(TestBoundingBoxInitialization)
{
  PRECICE_TEST();

  std::vector<Eigen::Vector3d> positions;
  std::vector<Eigen::Vector3d> data;
  std::vector<Eigen::Vector3d> expectedData;

  Eigen::Vector3d position;
  Eigen::Vector3d datum;

  for (int i = 0; i < 4; i++) {
    position[0] = i * 1.0;
    position[1] = i * 0.1;
    position[2] = -i * 10.0;
    positions.push_back(position);
    datum[0] = i * 1.0;
    datum[1] = i * 2.0;
    datum[2] = i * 3.0;
    data.push_back(datum);
    datum[0] = i * 1.0;
    datum[1] = i * 2.0;
    datum[2] = i * 3.0;
    expectedData.push_back(datum);
  }

  int i1 = -1, i2 = -1; // indices for data and positions

  if (context.isNamed("Fluid")) {
    if (context.isPrimary()) {
      i1 = 2;
      i2 = 4;
    } else {
      i1 = 0;
      i2 = 2;
    }
  } else {
    BOOST_TEST(context.isNamed("Structure"));
    // This participant starts with negated data
    for (int i = 0; i < 4; i++) {
      data[i] = -data[i];
    }
    if (context.isPrimary()) {
      i1 = 0;
      i2 = 2;
    } else {
      i1 = 2;
      i2 = 4;
    }
  }
  BOOST_REQUIRE(i1 >= 0);
  BOOST_REQUIRE(i2 >= 0);

  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  auto meshName = context.name + "Mesh";
  auto forcesID = "Forces";

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

  interface.advance(1.0);

  if (context.isNamed("Structure")) {
    double preciceDt = interface.getMaxTimeStepSize();
    for (size_t i = 0; i < vertexIDs.size(); i++) {
      interface.readData(meshName, forcesID, {&vertexIDs[i], 1}, preciceDt, {data[i + i1].data(), 3});
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
