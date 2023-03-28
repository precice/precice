#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(GlobalData)
// Test case with Global data for
// simple coupled simulation with iterations
// and without acceleration.
BOOST_AUTO_TEST_CASE(Implicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  // Set up Solverinterface
  precice::SolverInterface couplingInterface(context.name, context.config(), 0, 1);
  const int                dimensions = 2;
  BOOST_TEST(couplingInterface.getDimensions() == dimensions);

  std::string writeDataName;
  std::string readDataName;
  double      writeValue, expectedReadValue;

  if (context.isNamed("SolverOne")) {
    writeDataName     = "GlobalData1";
    readDataName      = "GlobalData2";
    writeValue        = 1;
    expectedReadValue = 2;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    writeDataName     = "GlobalData2";
    readDataName      = "GlobalData1";
    writeValue        = 2;
    expectedReadValue = 1;
  }
  int writeDataID = couplingInterface.getGlobalDataID(writeDataName);
  int readDataID  = couplingInterface.getGlobalDataID(readDataName);

  double              dt = 0;
  std::vector<double> writeData(dimensions, writeValue);
  std::vector<double> readData(dimensions, -1);

  // if (couplingInterface.requiresInitialData()) {
  //   BOOST_TEST(context.isNamed("SolverTwo"));
  //   couplingInterface.writeGlobalVectorData(writeDataID, writeData.data());
  // }

  dt = couplingInterface.initialize();

  while (couplingInterface.isCouplingOngoing()) {
    if (couplingInterface.requiresWritingCheckpoint()) {
    }
    // Write: from local data structure --> to precice buffer
    couplingInterface.writeGlobalVectorData(writeDataID, writeData.data());
    // Advance (exchange coupling data)
    dt = couplingInterface.advance(dt);
    // Read: from precice buffer --> to local data structure
    couplingInterface.readGlobalVectorData(readDataID, readData.data());
    // Check read data
    BOOST_TEST(expectedReadValue == readData.at(0));
    BOOST_TEST(expectedReadValue == readData.at(1));
    if (couplingInterface.requiresReadingCheckpoint()) {
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // GlobalData

#endif // PRECICE_NO_MPI
