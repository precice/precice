#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
// Test to ensure that RBF mappings can be used with 2LI
BOOST_AUTO_TEST_CASE(TestBoundingBoxInitializationRBF)
{
  PRECICE_TEST("Fluid"_on(2_ranks), "Structure"_on(2_ranks));

  std::vector<double> positions;
  std::vector<double> data;
  std::vector<double> expectedData;

  if (context.isNamed("Fluid")) {
    // Fluid 0 rank: create a grid from 0 to 20
    if (context.isPrimary()) {
      for (int x = 0; x < 20; ++x) {
        for (int y = 0; y < 2; ++y) {
          positions.push_back(x);
          positions.push_back(y);
          data.push_back(x + y * 7);
        }
      }
    } else {
      // Fluid 1 rank: create a grid from 20 to 40
      for (int x = 20; x < 40; ++x) {
        for (int y = 0; y < 2; ++y) {
          positions.push_back(x);
          positions.push_back(y);
          data.push_back(x + y * 7);
        }
      }
    }
  } else {
    BOOST_TEST(context.isNamed("Structure"));
    if (context.isPrimary()) {
      // Structure 0 rank: create a grid from 1 to 21
      for (int x = 1; x < 21; ++x) {
        for (int y = 0; y < 2; ++y) {
          positions.push_back(x);
          positions.push_back(y);
          // all vertices map to the 'last column'
          // of the primary fluid rank
          expectedData.push_back(12 + y * 7);
        }
      }
      data.resize(positions.size() / 2);
    } else {
      for (int x = 21; x < 41; ++x) {
        for (int y = 0; y < 2; ++y) {
          positions.push_back(x);
          positions.push_back(y);
          expectedData.push_back(x + y * 7);
        }
      }
      data.resize(positions.size() / 2);
    }
  }

  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  auto meshName = context.name + "Mesh";
  auto dataName = "Data";

  std::vector<int> vertexIDs(positions.size() / interface.getMeshDimensions(meshName));
  interface.setMeshVertices(meshName, positions, vertexIDs);

  interface.initialize();

  if (context.isNamed("Fluid")) {
    interface.writeData(meshName, dataName, vertexIDs, data);
  }

  interface.advance(1.0);

  if (context.isNamed("Structure")) {
    double preciceDt = interface.getMaxTimeStepSize();
    interface.readData(meshName, dataName, vertexIDs, preciceDt, data);
    // for (size_t i = 0; i < expectedData.size(); ++i) {
    //   BOOST_TEST(expectedData[i] == data[i]);
    // }
  }

  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI
