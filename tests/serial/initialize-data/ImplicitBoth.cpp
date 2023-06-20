#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(InitializeData)
/**
 * @brief Test simple coupled simulation with iterations, data initialization and without acceleration,
 * initialize=true for both participants. Tests https://github.com/precice/precice/issues/1367.
 *
 */
BOOST_AUTO_TEST_CASE(ImplicitBoth)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant couplingInterface(context.name, context.config(), 0, 1);

  std::string meshName;
  std::string writeDataName;
  std::string readDataName;
  double      writeValue, expectedReadValue;

  if (context.isNamed("SolverOne")) {
    meshName          = "MeshOne";
    writeDataName     = "Forces";
    readDataName      = "Velocities";
    writeValue        = 1;
    expectedReadValue = 2;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName          = "MeshTwo";
    writeDataName     = "Velocities";
    readDataName      = "Forces";
    writeValue        = 2;
    expectedReadValue = 1;
  }
  int                 dimensions = couplingInterface.getMeshDimensions(meshName);
  std::vector<double> vertex(dimensions, 0);
  int                 vertexID = couplingInterface.setMeshVertex(meshName, vertex);

  double              dt = 0;
  std::vector<double> writeData(dimensions, writeValue);
  std::vector<double> readData(dimensions, -1);

  if (couplingInterface.requiresInitialData()) {
    couplingInterface.writeData(meshName, writeDataName, {&vertexID, 1}, writeData);
  }

  couplingInterface.initialize();
  dt = couplingInterface.getMaxTimeStepSize();

  while (couplingInterface.isCouplingOngoing()) {
    if (couplingInterface.requiresWritingCheckpoint()) {
    }
    couplingInterface.readData(meshName, readDataName, {&vertexID, 1}, dt, readData);
    BOOST_TEST(expectedReadValue == readData.at(0));
    BOOST_TEST(expectedReadValue == readData.at(1));
    couplingInterface.writeData(meshName, writeDataName, {&vertexID, 1}, writeData);
    couplingInterface.advance(dt);
    dt = couplingInterface.getMaxTimeStepSize();
    if (couplingInterface.requiresReadingCheckpoint()) {
    }
  }
  couplingInterface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // InitializeData

#endif // PRECICE_NO_MPI
