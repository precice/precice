#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
// Test case for a direct mesh access on one participant to a mesh defined
// by another participant. The region of interest is defined through a
// boundingBox. As opposed to the 'boundingBoxExplicit' test case, this
// test case uses the same feature in an implicit setup.
BOOST_AUTO_TEST_CASE(Implicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  // Set up Solverinterface
  precice::SolverInterface couplingInterface(context.name, context.config(), 0, 1);
  BOOST_TEST(couplingInterface.getDimensions() == 2);
  constexpr int dim = 2;

  if (context.isNamed("SolverOne")) {
    std::vector<double>         positions   = {0.1, 0.1, 0.2, 0.05, 0.1, 0.0, 0.3, 0.9};
    std::array<double, dim * 2> boundingBox = {0.0, 1.0, 0.0, 1.0};
    std::vector<int>            ownIDs(4, -1);

    auto ownMeshName   = "MeshOne";
    auto otherMeshName = "MeshTwo";
    auto ownDataID     = "Forces";     //  ownMeshName
    auto otherDataID   = "Velocities"; //  otherMeshName

    // Define the own mesh
    couplingInterface.setMeshVertices(ownMeshName, ownIDs.size(), positions.data(), ownIDs.data());
    // TODO: Implement something in order to derive the bounding box from the mesh

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox.data());

    double dt = couplingInterface.initialize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int meshSize = couplingInterface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(meshSize == 3);

    // Allocate a vector containing the vertices
    std::vector<double> solverTwoMesh(meshSize * dim);
    std::vector<int>    otherIDs(meshSize);

    couplingInterface.getMeshVerticesAndIDs(otherMeshName, meshSize, otherIDs.data(), solverTwoMesh.data());
    // Some dummy writeData
    std::array<double, 3> writeData({1, 2, 3});

    // Expected data = positions of the other participant's mesh
    const std::vector<double> expectedData = {0.0, 0.0, 0.2, 0.3, 0.1, 0.1};
    BOOST_TEST(solverTwoMesh == expectedData);

    std::vector<double> readData(ownIDs.size(), -10);
    while (couplingInterface.isCouplingOngoing()) {
      if (couplingInterface.requiresWritingCheckpoint()) {
      }

      // Write data
      couplingInterface.writeBlockScalarData(otherMeshName, otherDataID, meshSize,
                                             otherIDs.data(), writeData.data());
      dt = couplingInterface.advance(dt);
      couplingInterface.readBlockScalarData(ownMeshName, ownDataID, ownIDs.size(),
                                            ownIDs.data(), readData.data());
      if (couplingInterface.requiresReadingCheckpoint()) {
      }

      // Expected data according to the writeData
      std::vector<double> expectedData({10, 11, 12, 13});
      BOOST_TEST(precice::testing::equals(expectedData, readData));
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    std::vector<double>         positions = {0.0, 0.0, 0.2, 0.3, 0.1, 0.1};
    std::vector<int>            ownIDs(3, -1);
    std::array<double, dim * 2> boundingBox = {0.0, 2.0, 0.0, 2.0};

    // Query IDs
    auto ownMeshName   = "MeshTwo";
    auto otherMeshName = "MeshOne";
    auto ownDataID     = "Velocities"; //  ownMeshName
    auto otherDataID   = "Forces";     //  otherMeshName

    // Define the mesh
    couplingInterface.setMeshVertices(ownMeshName, ownIDs.size(), positions.data(), ownIDs.data());
    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox.data());
    // Initialize
    double dt = couplingInterface.initialize();

    const int meshSize = couplingInterface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(meshSize == 4);

    // Allocate a vector containing the vertices
    std::vector<double> solverOneMesh(meshSize * dim);
    std::vector<int>    otherIDs(meshSize);

    couplingInterface.getMeshVerticesAndIDs(otherMeshName, meshSize, otherIDs.data(), solverOneMesh.data());
    // Some dummy writeData
    std::array<double, 4> writeData({10, 11, 12, 13});

    // Allocate data to read
    std::vector<double> readData(ownIDs.size(), -10);

    while (couplingInterface.isCouplingOngoing()) {
      if (couplingInterface.requiresWritingCheckpoint()) {
      }

      // Write data
      couplingInterface.writeBlockScalarData(otherMeshName, otherDataID, meshSize,
                                             otherIDs.data(), writeData.data());
      dt = couplingInterface.advance(dt);
      couplingInterface.readBlockScalarData(ownMeshName, ownDataID, ownIDs.size(),
                                            ownIDs.data(), readData.data());
      if (couplingInterface.requiresReadingCheckpoint()) {
      }

      // Expected data according to the writeData
      std::vector<double> expectedData({1, 2, 3});
      BOOST_TEST(precice::testing::equals(expectedData, readData));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // AccessReceivedMesh

#endif // PRECICE_NO_MPI
