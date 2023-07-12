#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(GlobalData)
// Test case to send global (meshless) data from one participant to another,
// while also writing data through direct mesh access.
// SolverOne writes all data and SolverTwo reads it.
BOOST_AUTO_TEST_CASE(ExplicitWithDirectMeshAccess)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  precice::Participant participant(context.name, context.config(), 0, 1);
  const int            dimensions = 2;

  // mesh related
  std::vector<double>                positions = {0.0, 0.0, 0.0, 0.05, 0.1, 0.1, 0.1, 0.0};
  std::vector<int>                   ids(4, -1);
  std::array<double, dimensions * 2> boundingBox = {0.0, 1.0, 0.0, 1.0};

  if (context.isNamed("SolverOne")) {
    auto otherMeshName = "MeshTwo";
    BOOST_REQUIRE(participant.getMeshDimensions(otherMeshName) == 2);
    const std::string globalScalarDataName = "GlobalScalarData";
    const std::string globalVectorDataName = "GlobalVectorData";
    const std::string dataName             = "ScalarData";

    // Define region of interest, where we could obtain direct write access
    participant.setMeshAccessRegion(otherMeshName, boundingBox);

    participant.initialize();
    double dt = participant.getMaxTimeStepSize();
    // Get the size of the filtered mesh within the bounding box
    // (provided by the coupling participant)
    const int meshSize = participant.getMeshVertexSize(otherMeshName);
    BOOST_TEST(meshSize == (ids.size()));
    // Allocate a vector containing the vertices
    std::vector<double> solverTwoMesh(meshSize * dimensions);
    participant.getMeshVerticesAndIDs(otherMeshName, ids, solverTwoMesh);
    // Expected data = positions of the other participant's mesh
    BOOST_TEST(precice::testing::equals(solverTwoMesh, positions));

    // Some dummy writeData
    std::vector<double> writeData{1, 2, 3, 4};
    double              writeGlobalScalarData{5};
    std::vector<double> writeGlobalVectorData(dimensions, 50.5);

    while (participant.isCouplingOngoing()) {
      // Write data to be sent to SolverTwo to buffer
      participant.writeData(otherMeshName, dataName, ids, writeData);
      participant.writeGlobalData(globalScalarDataName, {&writeGlobalScalarData, 1});
      participant.writeGlobalData(globalVectorDataName, writeGlobalVectorData);
      // send data
      participant.advance(dt);
      dt = participant.getMaxTimeStepSize();
      // change reference data for next check
      writeGlobalScalarData++;
      for (auto &elem : writeGlobalVectorData) {
        elem++;
      }
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));

    const std::string meshName             = "MeshTwo";
    const std::string dataName             = "ScalarData";
    const std::string globalScalarDataName = "GlobalScalarData";
    const std::string globalVectorDataName = "GlobalVectorData";
    BOOST_REQUIRE(participant.getMeshDimensions(meshName) == 2);

    // Define the mesh
    participant.setMeshVertices(meshName, positions, ids);

    // Allocate data to read
    std::vector<double> readData(4, std::numeric_limits<double>::max());
    double              readGlobalScalarData;
    std::vector<double> readGlobalVectorData(dimensions, -1);
    // Data expected to be received
    std::vector<double> expectedData{1, 2, 3, 4};
    double              expectedGlobalData{5};
    std::vector<double> expectedGlobalVectorData(dimensions, 50.5);

    // Initialize
    participant.initialize(); // For serial-explicit, first communication happens here
    double dt = participant.getMaxTimeStepSize();

    while (participant.isCouplingOngoing()) {
      // read received data from buffer
      participant.readData(meshName, dataName, ids, dt, readData);
      participant.readGlobalData(globalScalarDataName, dt, {&readGlobalScalarData, 1});
      participant.readGlobalData(globalVectorDataName, dt, readGlobalVectorData);
      // check if received data is correct
      BOOST_TEST(precice::testing::equals(expectedData, readData));
      BOOST_TEST(precice::testing::equals(expectedGlobalData, readGlobalScalarData));
      BOOST_TEST(precice::testing::equals(expectedGlobalVectorData, readGlobalVectorData));
      // receive next data
      participant.advance(dt);
      dt = participant.getMaxTimeStepSize();
      // change reference data for next check
      expectedGlobalData++;
      for (auto &elem : expectedGlobalVectorData) {
        elem++;
      }
      // Expected data according to the writeData
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // GlobalData

#endif // PRECICE_NO_MPI
