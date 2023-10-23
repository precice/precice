#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(TestBoundingBoxInitializationEmpty)
{
  PRECICE_TEST("Fluid"_on(2_ranks), "Structure"_on(2_ranks));

  std::vector<double> positions;
  std::vector<double> data;
  std::vector<double> expectedData;

  if (context.isNamed("Fluid")) {
    if (context.isPrimary()) {
      // init
    } else {
      // init
    }
  } else {
    BOOST_TEST(context.isNamed("Structure"));
    if (context.isPrimary()) {
      // init
    } else {
      // init
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
    // for (size_t d = 0; d < 3; d++) {
    //   BOOST_TEST(expectedData[i + i1][d] == data[i + i1][d]);
    // }
  }

  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI
