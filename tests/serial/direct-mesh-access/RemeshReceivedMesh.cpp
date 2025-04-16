#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test case for a direct mesh access by SolverTwo to a mesh defined
// by SolverOne. Both solvers read and write data to/from MeshOne.
// The received mesh is being reset
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(RemeshReceivedMesh)
{
  PRECICE_TEST();

  if (context.isNamed("SolverOne")) {
    // Set up Participant
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim              = 2;
    const auto           providedMeshName = "MeshOne";
    const auto           readDataName     = "Forces";
    const auto           writeDataName    = "Velocities";
    BOOST_TEST(interface.getMeshDimensions(providedMeshName) == 2);

    std::vector<double> positions = std::vector<double>({0.5, 0.25});
    int                 meshSize  = positions.size() / dim;
    std::vector<int>    ids(meshSize, -1);
    interface.setMeshVertices(providedMeshName, positions, ids);

    interface.initialize();

    // Some dummy writeData
    std::vector<double> readData(meshSize * dim, -1);
    std::vector<double> writeData(meshSize * dim, -1);
    std::vector<double> expectedData(meshSize * dim, 0);

    int time = 0;
    while (interface.isCouplingOngoing()) {

      time++;
      double dt = interface.getMaxTimeStepSize();

      // Read and compare against the reference
      interface.readData(providedMeshName, readDataName, ids, dt, readData);
      BOOST_TEST(expectedData == readData, boost::test_tools::per_element());

      // Now reset the mesh before we generate out write data
      // Not possible in the first time step, see #2093
      if (time > 1) {
        interface.resetMesh(providedMeshName);
        positions.push_back(0.2 + time * 0.1);
        positions.push_back(1.25);
        meshSize = positions.size() / dim;
        std::transform(positions.begin(), positions.end(), positions.begin(), [&](double value) {
          return value + 0.1;
        });

        ids.resize(meshSize, -1);
        interface.setMeshVertices(providedMeshName, positions, ids);
        writeData = readData = expectedData = std::vector<double>(meshSize * dim, -1);
      }

      // our artificial solve/ evaluation of the solve on our new mesh
      std::transform(positions.begin(), positions.end(), writeData.begin(), [&](double value) {
        return value * value - value * time;
      });
      interface.writeData(providedMeshName, writeDataName, ids, writeData);

      // triggers now the repartitioning and data exchange
      interface.advance(dt);

      // and the new reference data, according to the 'write' formula on the other participant (before we increment the time)
      std::transform(positions.begin(), positions.end(), expectedData.begin(), [&](double value) {
        return value * (value + 7.3) - (value + 5) * time;
      });
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    constexpr int        dim              = 2;
    const auto           receivedMeshName = "MeshOne";
    const auto           writeDataName    = "Forces";
    const auto           readDataName     = "Velocities";
    BOOST_TEST(interface.getMeshDimensions(receivedMeshName) == 2);

    std::array<double, dim * 2> boundingBox = std::array<double, dim * 2>{0.0, 1.0, 0.0, 2.0};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(receivedMeshName, boundingBox);

    interface.initialize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    int receivedMeshSize = interface.getMeshVertexSize(receivedMeshName);
    BOOST_TEST(receivedMeshSize == 1);

    // Allocate a vector containing the vertices
    std::vector<double> receivedMesh(receivedMeshSize * dim);
    std::vector<int>    receiveMeshIDs(receivedMeshSize, -1);
    interface.getMeshVertexIDsAndCoordinates(receivedMeshName, receiveMeshIDs, receivedMesh);

    // Allocate data to read and write
    std::vector<double> readData(receivedMeshSize * dim, -1);
    std::vector<double> writeData(receivedMeshSize * dim, -1);
    std::vector<double> expectedData(receivedMeshSize * dim, 0);

    // Expected data = positions of the other participant's mesh
    std::vector<double> expectedMesh({0.5, 0.25});
    BOOST_TEST(receivedMesh == expectedMesh, boost::test_tools::per_element());

    int time = 0;
    while (interface.isCouplingOngoing()) {

      time++;
      double dt = interface.getMaxTimeStepSize();

      // the reference data for this call will be generated based on the new coordinates
      // after we receive the new coordinates
      interface.readData(receivedMeshName, readDataName, receiveMeshIDs, dt, readData);
      BOOST_TEST(expectedData == readData, boost::test_tools::per_element());

      // our artificial solve
      std::transform(expectedMesh.begin(), expectedMesh.end(), writeData.begin(), [&](double value) {
        return value * (value + 7.3) - (value + 5) * time;
      });

      interface.writeData(receivedMeshName, writeDataName, receiveMeshIDs, writeData);
      interface.advance(dt);

      // mesh changes apply
      if (time > 1) {
        expectedMesh.push_back(0.2 + time * 0.1);
        expectedMesh.push_back(1.25);
        std::transform(expectedMesh.begin(), expectedMesh.end(), expectedMesh.begin(), [&](double value) {
          return value + 0.1;
        });
      }

      receivedMeshSize = interface.getMeshVertexSize(receivedMeshName);
      BOOST_TEST(receivedMeshSize == expectedMesh.size() / dim);
      receiveMeshIDs.resize(receivedMeshSize);
      writeData = readData = expectedData = receivedMesh = std::vector<double>(receivedMeshSize * dim, -1);
      interface.getMeshVertexIDsAndCoordinates(receivedMeshName, receiveMeshIDs, receivedMesh);
      BOOST_TEST(receivedMesh == expectedMesh, boost::test_tools::per_element());

      // and the new reference data, according to the 'write' formula on the other participant (before we increment the time)
      std::transform(expectedMesh.begin(), expectedMesh.end(), expectedData.begin(), [&](double value) {
        return value * value - value * time;
      });
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // DirectMeshAccess

#endif // PRECICE_NO_MPI
