#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
// Test case for a direct mesh access on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. The test case here is the most basic variant in order
// use such a feature. SolverTwo defines the mesh whereas SolverOne writes
// directly on this mesh.
BOOST_AUTO_TEST_CASE(Explicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  // Set up Solverinterface
  precice::SolverInterface couplingInterface(context.name, context.config(), 0, 1);
  BOOST_TEST(couplingInterface.getDimensions() == 2);

  std::vector<double> positions = {0.0, 0.0, 0.0, 0.05, 0.1, 0.1, 0.1, 0.0};
  std::vector<int>    ids(4, -1);

  constexpr int               dim         = 2;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0};

  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    auto dataName      = "Velocities"; //  otherMeshName

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox.data());

    double dt = couplingInterface.initialize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int meshSize = couplingInterface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(meshSize == (ids.size()));

    // Allocate a vector containing the vertices
    std::vector<double> solverTwoMesh(meshSize * dim);
    couplingInterface.getMeshVerticesAndIDs(otherMeshName, meshSize, ids.data(), solverTwoMesh.data());
    // Some dummy writeData
    std::array<double, 4> writeData({1, 2, 3, 4});

    // Expected data = positions of the other participant's mesh
    const std::vector<double> expectedData = positions;
    BOOST_TEST(precice::testing::equals(solverTwoMesh, expectedData));

    while (couplingInterface.isCouplingOngoing()) {
      // Write data
      couplingInterface.writeBlockScalarData(otherMeshName, dataName, meshSize,
                                             ids.data(), writeData.data());
      dt = couplingInterface.advance(dt);
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName = "MeshTwo";
    auto dataName = "Velocities"; //  meshName

    // Define the mesh
    couplingInterface.setMeshVertices(meshName, ids.size(), positions.data(), ids.data());
    // Allocate data to read
    std::vector<double> readData(4, std::numeric_limits<double>::max());

    // Initialize
    double dt = couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()) {

      dt = couplingInterface.advance(dt);
      couplingInterface.readBlockScalarData(meshName, dataName, ids.size(),
                                            ids.data(), readData.data());
      // Expected data according to the writeData
      std::vector<double> expectedData({1, 2, 3, 4});
      BOOST_TEST(precice::testing::equals(expectedData, readData));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // AccessReceivedMesh

#endif // PRECICE_NO_MPI
