#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(JustInTimeMapping)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

// Combining just-in-time mapping with direct mesh access
// One point is filtered out due to the defined mesh access region and still used for the jit mapping
// nearest-neighbor-consistent read and direct access write
BOOST_AUTO_TEST_CASE(MappingWithDirectAccessWrite)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 2;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, -0.011, 1.0};

  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    auto readDataName  = "Velocities"; // via jit
    auto writeDataName = "Mass-Flux";  // via direct-access
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 2);

    // Define region of interest, where we could obtain direct write access
    couplingInterface.setMeshAccessRegion(otherMeshName, boundingBox);

    couplingInterface.initialize();

    auto meshSize = couplingInterface.getMeshVertexSize(otherMeshName);
    BOOST_TEST(meshSize == 4);
    std::vector<double> receivedCoordinates(meshSize * dim);
    std::vector<int>    receivedIDs(meshSize);
    couplingInterface.getMeshVertexIDsAndCoordinates(otherMeshName, receivedIDs, receivedCoordinates);
    std::vector<double> expectedCoordinates({0.0, 0.0, 0.0, 0.05, 0.1, 0.1, 0.1, 0.0});
    BOOST_TEST(receivedCoordinates == expectedCoordinates, boost::test_tools::per_element());

    std::vector<double> writeData1({2.3, 2.4, 2.5, 2.6});
    std::vector<double> writeData2({23, 24, 25, 26});

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // read data:
      std::vector<double> expectedData1({5, 2, 3, 4, 2, 4});
      std::vector<double> expectedData2({-14, -11, -12, -13, -11, -13});
      std::vector<double> tmpPositions = {0.0, -0.01, 0.01, 0.05, 0.1, 0.1, 0.1, 0.0, 0, 0.05, 0.1, 0.0};

      for (std::size_t i = 0; i < expectedData1.size(); ++i) {
        std::vector<double> solverTwoCoord(dim);
        double              value;
        for (int d = 0; d < dim; ++d) {
          solverTwoCoord[d] = tmpPositions[i * dim + d];
        }

        couplingInterface.mapAndReadData(otherMeshName, readDataName, solverTwoCoord, dt, {&value, 1});
        // Expected data according to the writeData
        if (time == 1) {
          BOOST_TEST(expectedData1[i] == value);
        } else if (time == 2) {
          BOOST_TEST(expectedData2[i] == value);
        } else {
          PRECICE_ASSERT(false);
        }
      }

      // solve time step
      // write data (direct mesh access here)
      if (time == 1) {
        couplingInterface.writeData(otherMeshName, writeDataName, receivedIDs, writeData1);
      } else if (time == 2) {
        couplingInterface.writeData(otherMeshName, writeDataName, receivedIDs, writeData2);
      } else {
        PRECICE_ASSERT(false);
      }
      couplingInterface.advance(dt);
    }
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName      = "MeshTwo";
    auto writeDataName = "Velocities";
    auto readDataName  = "Mass-Flux";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName));

    // the last point is filtered out for direct access (not within the access region), but still considered for the jit mapping
    std::vector<double> positions = {0.0, 0.0, 0.0, 0.05, 0.1, 0.1, 0.1, 0.0, 0.0, -0.012};
    std::vector<int>    ids(5, -1);

    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // Some dummy readData
    std::array<double, 5> writeData1({1, 2, 3, 4, 5});
    std::array<double, 5> writeData2({-10, -11, -12, -13, -14});

    std::array<double, 5> readData;
    std::array<double, 5> readData0({0, 0, 0, 0});
    std::array<double, 5> readData1({2.3, 2.4, 2.5, 2.6, 0});
    std::array<double, 5> readData2({23, 24, 25, 26});
    // Initialize
    couplingInterface.initialize();
    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;
      // read data
      couplingInterface.readData(meshName, readDataName, ids, dt, readData);
      if (time == 1) {
        // the initial data from the other participant
        BOOST_TEST(readData == readData0, boost::test_tools::per_element());
      } else if (time == 2) {
        BOOST_TEST(readData == readData1, boost::test_tools::per_element());
      } else {
        PRECICE_ASSERT(false);
      }
      // solve time step
      // write data:
      if (time == 1) {
        couplingInterface.writeData(meshName, writeDataName, ids, writeData1);
      } else if (time == 2) {
        couplingInterface.writeData(meshName, writeDataName, ids, writeData2);
      } else {
        PRECICE_ASSERT(false);
      }
      couplingInterface.advance(dt);
    }
    couplingInterface.readData(meshName, readDataName, ids, 0, readData);
    BOOST_TEST(readData == readData2, boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
