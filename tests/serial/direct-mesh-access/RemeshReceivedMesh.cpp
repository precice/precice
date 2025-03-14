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
    double dt = interface.getMaxTimeStepSize();

    // Some dummy writeData
    std::vector<double> readData(ids.size(), -1);
    std::vector<double> writeData;
    for (int i = 0; i < meshSize; ++i)
      writeData.emplace_back(i + 5);

    int                 time = 0;
    std::vector<double> expectedData(1, 0);
    while (interface.isCouplingOngoing()) {

      interface.readData(providedMeshName, readDataName, ids, dt, readData);
      BOOST_TEST(expectedData == readData, boost::test_tools::per_element());
      expectedData[0] = 50;

      // Not possible in the first time step, see #2093
      if (time > 0) {
        interface.resetMesh(providedMeshName);
        positions.push_back(0.2 + time * 0.1);
        positions.push_back(1.25);
        meshSize = positions.size() / dim;
        ids.resize(meshSize, -1);
        interface.setMeshVertices(providedMeshName, positions, ids);
        writeData.push_back(writeData.size() + 5);
        readData.push_back(0);
        if (time > 1)
          expectedData.push_back(50 + expectedData.size());
        else
          expectedData.push_back(0);
      }

      interface.writeData(providedMeshName, writeDataName, ids, writeData);
      dt = interface.getMaxTimeStepSize();
      interface.advance(dt);
      time++;
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
    double dt = interface.getMaxTimeStepSize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    int receivedMeshSize = interface.getMeshVertexSize(receivedMeshName);
    BOOST_TEST(receivedMeshSize == 1);

    // Allocate a vector containing the vertices
    std::vector<double> receivedMesh(receivedMeshSize * dim);
    std::vector<int>    receiveMeshIDs(receivedMeshSize, -1);
    interface.getMeshVertexIDsAndCoordinates(receivedMeshName, receiveMeshIDs, receivedMesh);

    // Allocate data to read and write
    std::vector<double> readData(receiveMeshIDs.size(), -1);
    std::vector<double> writeData;
    for (int i = 0; i < receivedMeshSize; ++i)
      writeData.emplace_back(i + 50);

    // Expected data = positions of the other participant's mesh
    std::vector<double> expectedMesh({0.5, 0.25});
    BOOST_TEST(receivedMesh == expectedMesh, boost::test_tools::per_element());

    int                 time = 0;
    std::vector<double> expectedData(1, 0);
    while (interface.isCouplingOngoing()) {

      interface.readData(receivedMeshName, readDataName, receiveMeshIDs, dt, readData);
      BOOST_TEST(expectedData == readData, boost::test_tools::per_element());
      expectedData[0] = 5;
      interface.writeData(receivedMeshName, writeDataName, receiveMeshIDs, writeData);

      // Here, we would have the right place for resetting the access region
      // if (time > 0) {
      //   interface.resetMesh(receivedMeshName);
      //   interface.setMeshAccessRegion(receivedMeshName, boundingBox);
      // }

      interface.advance(dt);
      time++;

      if (time > 1) {
        receivedMeshSize = interface.getMeshVertexSize(receivedMeshName);
        BOOST_TEST(receivedMeshSize == (time));
        receiveMeshIDs.resize(receivedMeshSize);
        receivedMesh.resize(receivedMeshSize * dim, -1);
        interface.getMeshVertexIDsAndCoordinates(receivedMeshName, receiveMeshIDs, receivedMesh);
        expectedMesh.push_back(0.2 + (time - 1) * 0.1);
        expectedMesh.push_back(1.25);
        BOOST_TEST(receivedMesh == expectedMesh, boost::test_tools::per_element());
        writeData.emplace_back(writeData.size() + 50);
        readData.emplace_back(0);
        expectedData.push_back(5 + expectedData.size());
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // DirectMeshAccess

#endif // PRECICE_NO_MPI
