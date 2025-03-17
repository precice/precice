#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test case for a direct mesh access by SolverTwo to a mesh defined
// by SolverOne. Both solvers read and write data to/from MeshOne.
// The access region is being reset
// to read and write data, we encode the position vector of SolverOne
// and transform it to data to ensure that positions and data are matching
// after redefining the access regions
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ResetMeshAccessRegion)
{
  PRECICE_TEST();

  // We need to pick from this vector in both solvers
  std::vector<double> positions({0.5, 0.25, 1.0, 0.75, 1.5, 1.25});
  if (context.isNamed("SolverOne")) {
    // Set up Participant
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim              = 2;
    const auto           providedMeshName = "MeshOne";
    const auto           readDataName     = "Forces";
    const auto           writeDataName    = "Velocities";
    BOOST_TEST(interface.getMeshDimensions(providedMeshName) == 2);

    int              meshSize = positions.size() / dim;
    std::vector<int> ids(meshSize, -1);
    interface.setMeshVertices(providedMeshName, positions, ids);

    interface.initialize();

    // Some dummy writeData
    std::vector<double> readData(meshSize * dim, -1);
    std::vector<double> writeData(meshSize * dim, -1);

    int                 time = 0;
    std::vector<double> expectedData(meshSize * dim, 0);
    while (interface.isCouplingOngoing()) {

      time++;
      double dt = interface.getMaxTimeStepSize();

      // read
      interface.readData(providedMeshName, readDataName, ids, dt, readData);
      BOOST_TEST(expectedData == readData, boost::test_tools::per_element());

      // expected Data (lags one dt behind due to coupling scheme)
      std::transform(positions.begin(), positions.end(), expectedData.begin(), [&](double value) {
        return value * value - value * time;
      });

      // We manually account for values excluded due to the access region
      if (time < 3) {
        expectedData[2 * dim] = expectedData[2 * dim + 1] = 0;
      } else if (time < 4) {
        expectedData[0 * dim] = expectedData[0 * dim + 1] = 0;
      } else {
        expectedData[0 * dim] = expectedData[0 * dim + 1] = 0;
        expectedData[1 * dim] = expectedData[1 * dim + 1] = 0;
      }

      // The artificial solve
      std::transform(positions.begin(), positions.end(), writeData.begin(), [&](double value) {
        return value * value * time;
      });

      // write
      interface.writeData(providedMeshName, writeDataName, ids, writeData);
      interface.advance(dt);
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim              = 2;
    const auto           receivedMeshName = "MeshOne";
    const auto           writeDataName    = "Forces";
    const auto           readDataName     = "Velocities";
    BOOST_TEST(interface.getMeshDimensions(receivedMeshName) == 2);

    std::array<double, dim * 2> boundingBox{{0.0, 1.0, 0.0, 1.0}};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(receivedMeshName, boundingBox);

    interface.initialize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    int receivedMeshSize = interface.getMeshVertexSize(receivedMeshName);
    BOOST_TEST(receivedMeshSize == 2);

    // Allocate a vector containing the vertices
    std::vector<double> receivedMesh(receivedMeshSize * dim);
    std::vector<int>    receiveMeshIDs(receivedMeshSize, -1);
    interface.getMeshVertexIDsAndCoordinates(receivedMeshName, receiveMeshIDs, receivedMesh);

    // Allocate data to read and write
    std::vector<double> readData(receivedMeshSize * dim, -1);
    std::vector<double> writeData(receivedMeshSize * dim, -1);

    // Expected data = positions of the other participant's mesh
    auto                startIter = positions.begin();
    auto                endIter   = positions.begin() + 4;
    std::vector<double> expectedMesh(startIter, endIter);
    BOOST_TEST(receivedMesh == expectedMesh, boost::test_tools::per_element());

    int                 time = 0;
    std::vector<double> expectedData(receivedMeshSize * dim, 0);

    // The usual order for the precice time loop is read, solve, write advance
    // The usual order with dynamic meshes would be read, solve, reset, setVertices, write advance
    // However, the latter doesn't work for direct mesh access, since after redefining the mesh access region,
    // we don't have the new mesh we want to write on. Thus, the order has to be
    // read, solve, write, resetMeshAccessRegion, setMeshAccessRegion, advance, getMeshVerticesAndIDs
    while (interface.isCouplingOngoing()) {
      time++;
      double dt = interface.getMaxTimeStepSize();

      // read
      interface.readData(receivedMeshName, readDataName, receiveMeshIDs, dt, readData);
      BOOST_TEST(expectedData == readData, boost::test_tools::per_element());

      // (artificial) solve: fill writeData
      std::transform(expectedMesh.begin(), expectedMesh.end(), writeData.begin(), [&](double value) {
        return value * value - value * time;
      });

      // write
      interface.writeData(receivedMeshName, writeDataName, receiveMeshIDs, writeData);

      // we can't call this function in the first time step, see #2093
      if (time > 1 && time < 4) {
        interface.resetMeshAccessRegion(receivedMeshName);
        // we move the bounding box by 0.5 in each direction
        std::transform(boundingBox.begin(), boundingBox.end(), boundingBox.begin(),
                       [](double x) { return x + 0.5; });
        // For each movement here, we move the position pointer once forward
        std::advance(startIter, dim);
        endIter      = positions.end();
        expectedMesh = std::vector<double>(startIter, endIter);
        interface.setMeshAccessRegion(receivedMeshName, boundingBox);
      }
      interface.advance(dt);

      // After advance, we have the new vertices and IDs of the new access region
      // First, check the sizes
      receivedMeshSize = interface.getMeshVertexSize(receivedMeshName);
      BOOST_TEST(receivedMeshSize == (expectedMesh.size() / dim));

      // Resize vectors
      receiveMeshIDs.resize(receivedMeshSize, -1);
      receivedMesh = writeData = readData = expectedData = std::vector<double>(receivedMeshSize * dim, -1);
      interface.getMeshVertexIDsAndCoordinates(receivedMeshName, receiveMeshIDs, receivedMesh);

      BOOST_TEST(receivedMesh == expectedMesh, boost::test_tools::per_element());

      // the expected data lags one behind, so we generate (the same) reference solution as writeData for SolverOne, one dt later
      std::transform(expectedMesh.begin(), expectedMesh.end(), expectedData.begin(), [&](double value) {
        return value * value * time;
      });
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // DirectMeshAccess

#endif // PRECICE_NO_MPI
