#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
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

  SolverInterface couplingInterface(context.name, context.config(), 0, 1);

  int         dimensions = couplingInterface.getDimensions();
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
  std::vector<double> vertex(dimensions, 0);
  int                 vertexID = couplingInterface.setMeshVertex(meshName, vertex.data());

  double              dt = 0;
  std::vector<double> writeData(dimensions, writeValue);
  std::vector<double> readData(dimensions, -1);

  if (couplingInterface.requiresInitialData()) {
    couplingInterface.writeVectorData(meshName, writeDataName, vertexID, writeData.data());
  }

  dt = couplingInterface.initialize();

  while (couplingInterface.isCouplingOngoing()) {
    if (couplingInterface.requiresWritingCheckpoint()) {
    }
    couplingInterface.readVectorData(meshName, readDataName, vertexID, dt, readData.data());
    BOOST_TEST(expectedReadValue == readData.at(0));
    BOOST_TEST(expectedReadValue == readData.at(1));
    couplingInterface.writeVectorData(meshName, writeDataName, vertexID, writeData.data());
    dt = couplingInterface.advance(dt);
    if (couplingInterface.requiresReadingCheckpoint()) {
    }
  }
  couplingInterface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // InitializeData

#endif // PRECICE_NO_MPI
