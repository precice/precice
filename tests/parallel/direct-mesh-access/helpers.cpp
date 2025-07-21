#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include <numeric>
#include "precice/precice.hpp"
#include "testing/Testing.hpp"

void runTestMultipleBoundingBoxes2D(const TestContext &context)
{
  if (context.isNamed("SolverOne")) {
    // Set up Participant
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim           = 2;
    auto                 otherMeshName = "MeshTwo";
    auto                 writeDataName = "Velocities";
    BOOST_TEST(interface.getMeshDimensions(otherMeshName) == 2);

    std::array<double, dim * 2> boundingBox = context.isPrimary() ? std::array<double, dim * 2>{0.0, 1.0, 0.0, 2.1} : std::array<double, dim * 2>{0.0, 1.0, 3.5, 4.5};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(otherMeshName, boundingBox);

    // Second bounding box
    boundingBox = context.isPrimary() ? std::array<double, dim * 2>{0.0, 1.0, 1.9, 3.5} : std::array<double, dim * 2>{0.0, 1.0, 4.5, 5.5};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(otherMeshName, boundingBox);

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int otherMeshSize = interface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(otherMeshSize == 3);

    // Allocate a vector containing the vertices
    std::vector<double> solverTwoMesh(otherMeshSize * dim);
    std::vector<int>    otherIDs(otherMeshSize, -1);
    // Here, we don't receive vertex -0.5,-0.5, it's filtered out
    interface.getMeshVertexIDsAndCoordinates(otherMeshName, otherIDs, solverTwoMesh);
    // Expected data = positions of the other participant's mesh
    const std::vector<double> expectedData = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 2.0, 0.0, 3.5}) : std::vector<double>({0.0, 3.5, 0.0, 4.0, 0.0, 5.0});
    BOOST_TEST(solverTwoMesh == expectedData, boost::test_tools::per_element());

    // Some dummy writeData
    std::vector<double> writeData;
    for (int i = 0; i < otherMeshSize; ++i)
      writeData.emplace_back(i + 5 + (10 * context.isPrimary()));

    while (interface.isCouplingOngoing()) {
      // Write data
      interface.writeData(otherMeshName, writeDataName, otherIDs, writeData);
      interface.advance(dt);
      dt = interface.getMaxTimeStepSize();
      // reading data is not requires
    }
  } else {
    // Query IDs
    auto meshName     = "MeshTwo";
    auto readDataName = "Velocities";

    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    const int            dim = interface.getMeshDimensions(meshName);
    BOOST_TEST(context.isNamed("SolverTwo"));
    std::vector<double> positions = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 2.0, -0.5, -0.5}) : std::vector<double>({0.0, 3.5, 0.0, 4.0, 0.0, 5.0});
    std::vector<int>    ids(positions.size() / dim, -1);

    // Define the mesh
    interface.setMeshVertices(meshName, positions, ids);
    // Allocate data to read
    std::vector<double> readData(ids.size(), -1);

    // Initialize
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    while (interface.isCouplingOngoing()) {

      interface.advance(dt);
      dt = interface.getMaxTimeStepSize();

      interface.readData(meshName, readDataName, ids, dt, readData);
      // Expected data according to the writeData
      // Values are summed up
      std::vector<double> expectedData = context.isPrimary() ? std::vector<double>({15, 16, 0}) : std::vector<double>({22, 6, 7});
      BOOST_TEST(expectedData == readData, boost::test_tools::per_element());
    }
  }
}

void runTestMultipleBoundingBoxes3D(const TestContext &context)
{
  if (context.isNamed("SolverOne")) {
    // Set up Participant
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim           = 3;
    auto                 otherMeshName = "MeshTwo";
    auto                 writeDataName = "Velocities";
    auto                 readDataName  = "Forces";
    BOOST_TEST(interface.getMeshDimensions(otherMeshName) == 3);

    std::array<double, dim * 2> boundingBox = context.isPrimary() ? std::array<double, dim * 2>{0., 2.5, 0., 2.5, 0., 2.5} // 27
                                                                  : std::array<double, dim * 2>{10., 12., 0., 1., 0., 1.}; // 0
    interface.setMeshAccessRegion(otherMeshName, boundingBox);
    boundingBox = context.isPrimary() ? std::array<double, dim * 2>{1.5, 7.5, 0, 2.5, 0., 2.5} : std::array<double, dim * 2>{7.5, 10, 2.5, 5.0, 2.5, 5.0};
    interface.setMeshAccessRegion(otherMeshName, boundingBox);

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int otherMeshSize = interface.getMeshVertexSize(otherMeshName);
    if (context.isPrimary()) {
      BOOST_TEST(otherMeshSize == 72);
    } else {
      BOOST_TEST(otherMeshSize == 8);
    }
    // Allocate a vector containing the vertices
    std::vector<double> solverTwoMesh(otherMeshSize * dim);
    std::vector<int>    otherIDs(otherMeshSize, -1);
    // Here, we don't receive vertex -0.5,-0.5, it's filtered out
    interface.getMeshVertexIDsAndCoordinates(otherMeshName, otherIDs, solverTwoMesh);
    // Expected data = positions of the other participant's mesh
    std::vector<double> expectedData;
    if (context.isPrimary()) {
      expectedData = std::vector<double>{
          /* z = 0 (y = 0‥5, x = 0‥5) */
          0, 0, 0, 1, 0, 0,
          2, 0, 0, 3, 0, 0,
          4, 0, 0, 0, 1, 0,
          1, 1, 0, 2, 1, 0,
          3, 1, 0, 4, 1, 0,
          0, 2, 0, 1, 2, 0,
          /* z = 1 */
          2, 2, 0, 3, 2, 0,
          4, 2, 0, 0, 0, 1,
          1, 0, 1, 2, 0, 1,
          3, 0, 1, 4, 0, 1,
          0, 1, 1, 1, 1, 1,
          2, 1, 1, 3, 1, 1,
          /* z = 2 */
          4, 1, 1, 0, 2, 1,
          1, 2, 1, 2, 2, 1,
          3, 2, 1, 4, 2, 1,
          0, 0, 2, 1, 0, 2,
          2, 0, 2, 3, 0, 2,
          4, 0, 2, 0, 1, 2,
          /* z = 3 */
          1, 1, 2, 2, 1, 2,
          3, 1, 2, 4, 1, 2,
          0, 2, 2, 1, 2, 2,
          2, 2, 2, 3, 2, 2,
          4, 2, 2, 5, 0, 0,
          6, 0, 0, 7, 0, 0,
          /* z = 4 */
          5, 1, 0, 6, 1, 0,
          7, 1, 0, 5, 2, 0,
          6, 2, 0, 7, 2, 0,
          5, 0, 1, 6, 0, 1,
          7, 0, 1, 5, 1, 1,
          6, 1, 1, 7, 1, 1,
          /* z = 5 */
          5, 2, 1, 6, 2, 1,
          7, 2, 1, 5, 0, 2,
          6, 0, 2, 7, 0, 2,
          5, 1, 2, 6, 1, 2,
          7, 1, 2, 5, 2, 2,
          6, 2, 2, 7, 2, 2};
    } else {
      expectedData = std::vector<double>{
          8, 3, 3, 9, 3, 3, 8, 4, 3, 9, 4, 3, 8, 3, 4, 9, 3, 4, 8, 4, 4, 9, 4, 4};
    }
    BOOST_TEST(solverTwoMesh == expectedData, boost::test_tools::per_element());

    // Some dummy writeData
    std::vector<double> writeData;
    for (int i = 0; i < otherMeshSize; ++i)
      writeData.emplace_back(i + 5 + (100 * context.isPrimary()));

    while (interface.isCouplingOngoing()) {
      // Write data
      interface.writeData(otherMeshName, writeDataName, otherIDs, writeData);
      interface.advance(dt);
      dt = interface.getMaxTimeStepSize();
      // reading data is not requires
    }
  } else {
    // Query IDs
    auto meshName      = "MeshTwo";
    auto readDataName  = "Velocities";
    auto writeDataName = "Forces";

    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    const int            dim = interface.getMeshDimensions(meshName);
    BOOST_TEST(context.isNamed("SolverTwo"));
    std::vector<double> positions;

    if (context.isPrimary()) {
      const int nodesPerEdge = 5;
      for (int z = 0; z < nodesPerEdge; ++z)
        for (int y = 0; y < nodesPerEdge; ++y)
          for (int x = 0; x < nodesPerEdge; ++x) {
            positions.push_back({static_cast<double>(x)});
            positions.push_back({static_cast<double>(y)});
            positions.push_back({static_cast<double>(z)});
          }
    } else {
      const int nodesPerEdge = 5;
      for (int z = 0; z < nodesPerEdge; ++z)
        for (int y = 0; y < nodesPerEdge; ++y)
          for (int x = 5; x < nodesPerEdge + 5; ++x) { // shifted by 5
            positions.push_back({static_cast<double>(x)});
            positions.push_back({static_cast<double>(y)});
            positions.push_back({static_cast<double>(z)});
          }
    }
    std::vector<int> ids(positions.size() / dim, -1);
    // Define the mesh
    interface.setMeshVertices(meshName, positions, ids);
    // Allocate data to read
    std::vector<double> readData(ids.size(), -1);

    // Initialize
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    while (interface.isCouplingOngoing()) {

      interface.advance(dt);
      dt = interface.getMaxTimeStepSize();

      interface.readData(meshName, readDataName, ids, dt, readData);
    }
  }
}

// StartIndex is here the first index to be used for writing on the secondary rank
void runTestAccessReceivedMesh(const TestContext        &context,
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
                              (!context.isPrimary() && meshSize == static_cast<int>(expectedPositionSecondaryRank.size()) / dim);
    BOOST_TEST(expectedSize);

    // Allocate memory
    std::vector<int>    ids(meshSize);
    std::vector<double> coordinates(static_cast<std::size_t>(meshSize) * static_cast<std::size_t>(dim));
    interface.getMeshVertexIDsAndCoordinates(otherMeshName, ids, coordinates);

    // Check the received vertex coordinates
    std::vector<double> expectedPositions = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 2.0, 0.0, 3.0}) : expectedPositionSecondaryRank;
    BOOST_TEST(testing::equals(expectedPositions, coordinates));

    // Check the received vertex IDs (IDs are local?!)
    std::vector<int> expectedIDs;
    for (std::size_t i = 0; i < static_cast<size_t>(meshSize); ++i)
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
        // Corresponds semantically to meshSize - startIndex > 0
        // but meshSize - startIndex > 0 might underflow and the static analysis complained
        if (meshSize > static_cast<int>(startIndex)) {
          const int *ids_ptr  = &ids.at(startIndex);
          const auto vertices = meshSize - startIndex;
          interface.writeData(otherMeshName, dataName, {ids_ptr, vertices}, {writeData.data(), vertices});
        }
      }

      interface.advance(dt);
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

    // This is not allowed, as meshName is a local mesh
    std::vector<double> dummyBB({0.0, 1.0, 0.0, 3.5});
    BOOST_CHECK_THROW(interface.setMeshAccessRegion(meshName, dummyBB), ::precice::Error);

    {
      // Check, that we can't use the 'getMeshVertexIDsAndCoordinates' function on the provided meshes,
      // (the actual purpose is of course using it on received meshes)
      const std::size_t ownMeshSize = interface.getMeshVertexSize(meshName);
      BOOST_TEST(ownMeshSize == size);
      std::vector<int>    ownIDs(ownMeshSize);
      std::vector<double> ownCoordinates(ownMeshSize * dim);
      BOOST_CHECK_THROW(interface.getMeshVertexIDsAndCoordinates(meshName, ownIDs, ownCoordinates), ::precice::Error);
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
