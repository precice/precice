#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include <numeric>
#include "precice/precice.hpp"
#include "testing/Testing.hpp"

// StartIndex is here the first index to be used for writing on the secondary rank
void runTestAccessReceivedMesh(const TestContext &       context,
                               const std::vector<double> boundingBoxSecondaryRank,
                               const std::vector<double> writeDataSecondaryRank,
                               const std::vector<double> expectedPositionSecondaryRank,
                               const std::vector<double> expectedReadDataSecondaryRank,
                               const size_t              startIndex)
{
  if (context.isNamed("SolverOne")) {
    // Defines the bounding box and writes data to the received mesh
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    auto                 otherMeshName = "MeshTwo";
    auto                 dataName      = "Velocities";
    const int            dim           = interface.getMeshDimensions(otherMeshName);

    std::vector<double> boundingBox = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 3.5}) : boundingBoxSecondaryRank;
    // Set bounding box
    interface.setMeshAccessRegion(otherMeshName, boundingBox);
    // Initialize the Participant
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Get relevant size, allocate data structures and retrieve coordinates
    const auto meshSize = interface.getMeshVertexSize(otherMeshName);

    // According to the bounding boxes and vertices: the primary rank receives 3 vertices, the secondary rank 2
    const bool expectedSize = (context.isPrimary() && meshSize == 3) ||
                              (!context.isPrimary() && meshSize == expectedPositionSecondaryRank.size() / dim);
    BOOST_TEST(expectedSize);

    // Allocate memory
    std::vector<int>    ids(meshSize);
    std::vector<double> coordinates(meshSize * dim);
    interface.getMeshVertexIDsAndCoordinates(otherMeshName, ids, coordinates);

    // Check the received vertex coordinates
    std::vector<double> expectedPositions = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 2.0, 0.0, 3.0}) : expectedPositionSecondaryRank;
    BOOST_TEST(testing::equals(expectedPositions, coordinates));

    // Check the received vertex IDs (IDs are local?!)
    std::vector<int> expectedIDs;
    for (std::size_t i = 0; i < meshSize; ++i)
      expectedIDs.emplace_back(i);
    BOOST_TEST(expectedIDs == ids);

    // Create some unique writeData in order to check it in the other participant
    std::vector<double> primaryData(meshSize);
    std::iota(primaryData.begin(), primaryData.end(), 1);
    std::vector<double> writeData = context.isPrimary() ? primaryData : writeDataSecondaryRank;

    while (interface.isCouplingOngoing()) {
      // Write data
      if (context.isPrimary()) {
        interface.writeData(otherMeshName, dataName, ids, writeData);
      } else {
        if (meshSize - startIndex > 0) {
          const int *ids_ptr  = &ids.at(startIndex);
          const auto vertices = meshSize - startIndex;
          interface.writeData(otherMeshName, dataName, {ids_ptr, vertices}, {writeData.data(), vertices});
        }
      }

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
    }
  } else {
    // Defines the mesh and reads data
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    precice::Participant interface(context.name, context.config(), context.rank, context.size);

    // Get IDs
    auto      meshName = "MeshTwo";
    auto      dataName = "Velocities";
    const int dim      = interface.getMeshDimensions(meshName);
    BOOST_TEST(dim == 2);
    // Define the interface
    std::vector<double> positions = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 2.0}) : std::vector<double>({0.0, 3.0, 0.0, 4.0, 0.0, 5.0});

    const int        size = positions.size() / dim;
    std::vector<int> ids(size);

    interface.setMeshVertices(meshName, positions, ids);

    {
      // Check, if we can use the 'getMeshVertexIDsAndCoordinates' function on provided meshes as well,
      // though the actual purpose is of course using it on received meshes
      const std::size_t ownMeshSize = interface.getMeshVertexSize(meshName);
      BOOST_TEST(ownMeshSize == size);
      std::vector<int>    ownIDs(ownMeshSize);
      std::vector<double> ownCoordinates(ownMeshSize * dim);
      interface.getMeshVertexIDsAndCoordinates(meshName, ownIDs, ownCoordinates);
      BOOST_TEST(ownIDs == ids);
      BOOST_TEST(testing::equals(positions, ownCoordinates));
    }

    // Initialize the Participant
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Start the time loop
    std::vector<double> readData(size);
    while (interface.isCouplingOngoing()) {

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
      interface.readData(meshName, dataName, ids, dt, readData);

      // Check the received data
      const std::vector<double> expectedReadData = context.isPrimary() ? std::vector<double>({1, 2}) : expectedReadDataSecondaryRank;
      BOOST_TEST(expectedReadData == readData, boost::test_tools::per_element());
    }
  }
}

#endif
