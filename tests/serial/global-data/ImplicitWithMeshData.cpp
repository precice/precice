#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(GlobalData)
// Test case with Global data for
// simple coupled simulation with iterations
// and without acceleration.
BOOST_AUTO_TEST_CASE(ImplicitWithMeshData)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  // Set up API
  precice::Participant participant(context.name, context.config(), 0, 1);

  std::string meshName;
  std::string writeDataName;
  std::string readDataName;
  std::string writeGlobalDataName;
  std::string readGlobalDataName;
  double      writeValue, expectedReadValue, writeGlobalValue, expectedReadGlobalValue;

  if (context.isNamed("SolverOne")) {
    // mesh data
    meshName          = "MeshOne";
    writeDataName     = "Forces";
    readDataName      = "Velocities";
    writeValue        = 1;
    expectedReadValue = 2;
    // global data
    writeGlobalDataName     = "GlobalData1";
    readGlobalDataName      = "GlobalData2";
    writeGlobalValue        = 3;
    expectedReadGlobalValue = 4;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // mesh data
    meshName          = "MeshTwo";
    writeDataName     = "Velocities";
    readDataName      = "Forces";
    writeValue        = 2;
    expectedReadValue = 1;
    // global data
    writeGlobalDataName     = "GlobalData2";
    readGlobalDataName      = "GlobalData1";
    writeGlobalValue        = 4;
    expectedReadGlobalValue = 3;
  }

  // create mesh with one vertex
  const int           dimensions = participant.getMeshDimensions(meshName);
  std::vector<double> vertex(dimensions, 0);
  int                 vertexID = participant.setMeshVertex(meshName, vertex);

  double              dt = 0;
  std::vector<double> writeData(dimensions, writeValue);
  std::vector<double> readData(dimensions, -1);
  std::vector<double> writeGlobalData(dimensions, writeGlobalValue);
  std::vector<double> readGlobalData(dimensions, -1);

  if (participant.requiresInitialData()) {
    BOOST_TEST(context.isNamed("SolverTwo"));
    participant.writeData(meshName, writeDataName, {&vertexID, 1}, writeData);
  }

  participant.initialize();
  dt = participant.getMaxTimeStepSize();

  while (participant.isCouplingOngoing()) {
    if (participant.requiresWritingCheckpoint()) {
    }
    // Write global data
    participant.writeGlobalData(writeGlobalDataName, writeGlobalData);
    // Read mesh data
    participant.readData(meshName, readDataName, {&vertexID, 1}, dt, readData);
    BOOST_TEST(expectedReadValue == readData.at(0));
    BOOST_TEST(expectedReadValue == readData.at(1));
    // Write mesh data
    participant.writeData(meshName, writeDataName, {&vertexID, 1}, writeData);
    // Advance (exchange coupling data)
    participant.advance(dt);
    dt = participant.getMaxTimeStepSize();
    // Read global data
    participant.readGlobalData(readGlobalDataName, dt, readGlobalData);
    BOOST_TEST(expectedReadGlobalValue == readGlobalData.at(0));
    BOOST_TEST(expectedReadGlobalValue == readGlobalData.at(1));
    if (participant.requiresReadingCheckpoint()) {
    }
  }
  participant.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // GlobalData

#endif // PRECICE_NO_MPI
