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
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(DirectMeshAccess)
PRECICE_TEST_SETUP("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks))
BOOST_AUTO_TEST_CASE(ResetMeshAccessRegion)
{
  PRECICE_TEST();

  // We need to pick from this vector in both solvers
  constexpr int       dim = 2;
  std::vector<double> globalPositions(40 * dim);
  // Generate a vector from 0 to 40
  std::generate(globalPositions.begin(), globalPositions.end(), [n = 1, i = 0, &context]() mutable { return (i++ % dim == 0 ? n++ : 0); });

  if (context.isNamed("SolverOne")) {
    // for rank zero x coords 1..20 and rank 1 x coords 21...40
    precice::span<double> positions(globalPositions.data() + 20 * dim * context.rank, 20 * dim);
    // Set up Participant
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
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
        return value * value - value * (time);
      });

      // There is a duplicate write in the overlap region
      for (std::size_t i = 0; i < meshSize; ++i) {
        auto start = time > 1 ? 10 : 10;
        auto end   = time > 1 ? 31 : 31;
        if (positions[i * dim] > start &&
            positions[i * dim] < end) {
          expectedData[i * dim] *= 2;
        }
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

    std::array<double, dim * 2> boundingBox{{0.5 + context.rank * 10, 30.5 + context.rank * 20, 0.0, 1.0}};
    // Define region of interest, where we could obtain direct write access
    interface.setMeshAccessRegion(receivedMeshName, boundingBox);

    interface.initialize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    int receivedMeshSize = interface.getMeshVertexSize(receivedMeshName);
    BOOST_TEST(receivedMeshSize == 30);

    // Allocate a vector containing the vertices
    std::vector<double> receivedMesh(receivedMeshSize * dim);
    std::vector<int>    receiveMeshIDs(receivedMeshSize, -1);
    interface.getMeshVertexIDsAndCoordinates(receivedMeshName, receiveMeshIDs, receivedMesh);

    // Allocate data to read and write
    std::vector<double> readData(receivedMeshSize * dim, -1);
    std::vector<double> writeData(receivedMeshSize * dim, -1);

    // Expected data = positions of the other participant's mesh
    auto                startIter = context.rank == 0 ? globalPositions.begin() : globalPositions.begin() + context.rank * 10 * dim;
    auto                endIter   = context.rank == 0 ? startIter + dim * 30 : globalPositions.end();
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

      if (context.rank == 0) {
        for (int i = 0; i < expectedMesh.size() / dim; ++i)
          std::cout << "Writing: Rank 0: " << writeData[dim * i] << "  positions: " << expectedMesh[dim * i] << "  at time: " << time << std::endl;
      }
      if (context.rank == 1) {
        for (int i = 0; i < expectedMesh.size() / dim; ++i)
          std::cout << "Writing: Rank 1: " << writeData[dim * i] << "  positions: " << expectedMesh[dim * i] << "  at time: " << time << std::endl;
      }

      // write
      interface.writeData(receivedMeshName, writeDataName, receiveMeshIDs, writeData);

      // we can't call this function in the first time step, see #2093
      if (time > 1) {
        interface.resetMeshAccessRegion(receivedMeshName);
        // we move the bounding box in the middle region
        // std::transform(boundingBox.begin(), boundingBox.begin() + 2, boundingBox.begin(),
        //                [](double x) { return x + 1; });
        // For each movement here, we move the position pointer once forward
        if (context.rank == 0) {
          std::advance(endIter, dim);
          boundingBox[1] += 1;
        }
        if (context.rank == 1) {
          std::advance(startIter, dim);
          boundingBox[0] += 1;
        }
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
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // DirectMeshAccess

#endif // PRECICE_NO_MPI
