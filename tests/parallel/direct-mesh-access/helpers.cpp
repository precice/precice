#include <cstddef>
#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "precice/SolverInterface.hpp"
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
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    auto                     otherMeshID = "MeshTwo";
    auto                     dataID      = "Velocities"; //  otherMeshID
    const int                dim         = interface.getDimensions();

    std::vector<double> boundingBox = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 3.5}) : boundingBoxSecondaryRank;
    // Set bounding box
    interface.setMeshAccessRegion(otherMeshID, boundingBox.data());
    // Initialize the solverinterface
    double dt = interface.initialize();

    // Get relevant size, allocate data structures and retrieve coordinates
    const std::size_t meshSize = interface.getMeshVertexSize(dataID);

    // According to the bounding boxes and vertices: the primary rank receives 3 vertices, the secondary rank 2
    const bool expectedSize = (context.isPrimary() && meshSize == 3) ||
                              (!context.isPrimary() && meshSize == expectedPositionSecondaryRank.size() / dim);
    BOOST_TEST(expectedSize);

    // Allocate memory
    std::vector<int>    ids(meshSize);
    std::vector<double> coordinates(meshSize * dim);
    interface.getMeshVerticesAndIDs(otherMeshID, meshSize, ids.data(), coordinates.data());

    // Check the received vertex coordinates
    std::vector<double> expectedPositions = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 2.0, 0.0, 3.0}) : expectedPositionSecondaryRank;
    BOOST_TEST(testing::equals(expectedPositions, coordinates));

    // Check the received vertex IDs (IDs are local?!)
    std::vector<int> expectedIDs;
    for (std::size_t i = 0; i < meshSize; ++i)
      expectedIDs.emplace_back(i);
    BOOST_TEST(expectedIDs == ids);

    // Create some unique writeData in order to check it in the other participant
    std::vector<double> writeData = context.isPrimary() ? std::vector<double>({1, 2, 3}) : writeDataSecondaryRank;

    while (interface.isCouplingOngoing()) {
      // Write data
      if (context.isPrimary()) {
        interface.writeBlockScalarData(meshID, dataID, meshSize,
                                       ids.data(), writeData.data());
      } else {
        // In order to prevent hypothetical index overruns reported by glibcc
        const int *ids_ptr = startIndex < ids.size() ? &ids[startIndex] : nullptr;
        interface.writeBlockScalarData(meshID, dataID, meshSize - startIndex,
                                       ids_ptr, writeData.data());
      }

      dt = interface.advance(dt);
    }
  } else {
    // Defines the mesh and reads data
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    BOOST_TEST(interface.getDimensions() == 2);

    // Get IDs
    auto      meshID = "MeshTwo";
    auto      dataID = "Velocities"; //  meshID
    const int dim    = interface.getDimensions();
    // Define the interface
    std::vector<double> positions = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 2.0}) : std::vector<double>({0.0, 3.0, 0.0, 4.0, 0.0, 5.0});

    const int        size = positions.size() / dim;
    std::vector<int> ids(size);

    interface.setMeshVertices(meshID, size, positions.data(), ids.data());

    {
      // Check, if we can use the 'getMeshVerticesAndIDs' function on provided meshes as well,
      // though the actual purpose is of course using it on received meshes
      const std::size_t ownMeshSize = interface.getMeshVertexSize(meshID);
      BOOST_TEST(ownMeshSize == size);
      std::vector<int>    ownIDs(ownMeshSize);
      std::vector<double> ownCoordinates(ownMeshSize * dim);
      interface.getMeshVerticesAndIDs(meshID, ownMeshSize, ownIDs.data(), ownCoordinates.data());
      BOOST_TEST(ownIDs == ids);
      BOOST_TEST(testing::equals(positions, ownCoordinates));
    }

    // Initialize the solverinterface
    double dt = interface.initialize();

    // Start the time loop
    std::vector<double> readData(size);
    while (interface.isCouplingOngoing()) {

      dt = interface.advance(dt);
      interface.readBlockScalarData(meshID, dataID, size,
                                    ids.data(), readData.data());

      // Check the received data
      const std::vector<double> expectedReadData = context.isPrimary() ? std::vector<double>({1, 2}) : expectedReadDataSecondaryRank;
      BOOST_TEST(expectedReadData == readData);
    }
  }
}

#endif
