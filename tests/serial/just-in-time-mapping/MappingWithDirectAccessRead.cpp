#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(JustInTimeMapping)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

// Combining just-in-time mapping with direct mesh access
// nearest-neighbor-conservative write and direct access read
BOOST_AUTO_TEST_CASE(MappingWithDirectAccessRead)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant couplingInterface(context.name, context.config(), 0, 1);

  constexpr int               dim         = 2;
  std::array<double, dim * 2> boundingBox = {0.0, 1.0, -0.011, 1.0};

  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    auto writeDataName = "Velocities";
    auto readDataName  = "Mass-Flow";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(otherMeshName) == 2);

    std::vector<double> writeData1({1, 2, -3, 4, 2});
    std::vector<double> writeData2({-4, 12, 3, 5, 7});

    std::array<double, 4> expectedData0({0, 0, 0, 0});
    std::array<double, 4> expectedData1({1, 2, 3, 4});
    std::array<double, 4> expectedData2({-10, -11, -12, -13});
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
    std::vector<double> readData(meshSize);

    double time = 0;
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;

      // read data
      couplingInterface.readData(otherMeshName, readDataName, receivedIDs, dt, readData);
      if (time == 1) {
        BOOST_TEST(readData == expectedData0, boost::test_tools::per_element());
      } else if (time == 2) {
        BOOST_TEST(readData == expectedData1, boost::test_tools::per_element());
      } else {
        BOOST_TEST(false);
      }

      // solve time step
      // write data:
      std::vector<double> tmpPositions = {0.0, -0.01, 0.01, 0.05, 0.1, 0.1, 0.1, 0.0, 0, 0.05};
      for (std::size_t i = 0; i < writeData1.size(); ++i) {
        std::vector<double> solverTwoCoord(dim);
        for (int d = 0; d < dim; ++d) {
          solverTwoCoord[d] = tmpPositions[i * dim + d];
        }
        if (time == 1) {
          couplingInterface.writeAndMapData(otherMeshName, writeDataName, solverTwoCoord, {&writeData1[i], 1});
        } else if (time == 2) {
          couplingInterface.writeAndMapData(otherMeshName, writeDataName, solverTwoCoord, {&writeData2[i], 1});
        } else {
          BOOST_TEST(false);
        }
      }
      couplingInterface.advance(dt);
    }
    couplingInterface.readData(otherMeshName, readDataName, receivedIDs, 0, readData);
    BOOST_TEST(readData == expectedData2, boost::test_tools::per_element());
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    auto meshName      = "MeshTwo";
    auto readDataName  = "Velocities";
    auto writeDataName = "Mass-Flow";
    BOOST_REQUIRE(couplingInterface.getMeshDimensions(meshName));

    std::vector<double> positions = {0.0, 0.0, 0.0, 0.05, 0.1, 0.1, 0.1, 0.0, 0.0, -0.012};
    std::vector<int>    ids(5, -1);

    // the last value maps to a coordinate not within the bounding box of the jit mapping
    std::vector<double> expectedData1({0, 4, -3, 4, 1});
    std::vector<double> expectedData2({0, 19, 3, 5, -4});

    // Define the mesh
    couplingInterface.setMeshVertices(meshName, positions, ids);
    // Some dummy readData
    std::vector<double> readData(5, -1);

    std::array<double, 5> writeData1({1, 2, 3, 4, 5});
    std::array<double, 5> writeData2({-10, -11, -12, -13, -14});

    double time = 0;
    // Initialize
    couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()) {
      double dt = couplingInterface.getMaxTimeStepSize();
      time += dt;
      couplingInterface.readData(meshName, readDataName, ids, dt, readData);

      if (time == 1) {
        BOOST_TEST(readData == expectedData1, boost::test_tools::per_element());
      } else if (time == 2) {
        BOOST_TEST(readData == expectedData2, boost::test_tools::per_element());
      } else {
        BOOST_TEST(false);
      }
      // solve time step
      // write data
      if (time == 1) {
        couplingInterface.writeData(meshName, writeDataName, ids, writeData1);
      } else if (time == 2) {
        couplingInterface.writeData(meshName, writeDataName, ids, writeData2);
      } else {
        BOOST_TEST(false);
      }
      couplingInterface.advance(dt);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
