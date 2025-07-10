#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test case for a direct mesh access by SolverTwo to a mesh defined
// by SolverOne.
// Here, SolverOne defines positions {0.5, 0.25, 1.0, 0.75, 1.5, 1.25}
// SolverTwo defines the access region {0.5, 1.5, 0.5, 1.5}, which filters
// out the first vertex, reading and writing (should) just happen on the
// latter two
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(DirectAccessFiltered)
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
      // BOOST_TEST(expectedData == readData, boost::test_tools::per_element());

      // expected Data (lags one dt behind due to coupling scheme)
      std::transform(positions.begin(), positions.end(), expectedData.begin(), [&](double value) {
        return value * value - value * time;
      });

      // We manually account for values excluded due to the access region
      expectedData[2 * dim] = expectedData[2 * dim + 1] = 0;

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

    // This bounding box excludes the first vertex of the other participant
    std::array<double, dim * 2> boundingBox{{0.5, 1.5, 0.5, 1.5}};
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
    auto                startIter = positions.begin() + 2;
    auto                endIter   = positions.begin() + 6;
    std::vector<double> expectedMesh(startIter, endIter);
    BOOST_TEST(receivedMesh == expectedMesh, boost::test_tools::per_element());

    int                 time = 0;
    std::vector<double> expectedData(receivedMeshSize * dim, 0);

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

      interface.writeData(receivedMeshName, writeDataName, receiveMeshIDs, writeData);
      interface.advance(dt);

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
