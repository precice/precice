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

  std::string writeVectorDataName, writeScalarDataName;
  std::string readVectorDataName, readScalarDataName;
  double      writeVectorValue, expectedReadVectorValue, writeScalarValue, expectedReadScalarValue;

  if (context.isNamed("SolverOne")) {
    writeVectorDataName     = "GlobalVectorData1";
    writeScalarDataName     = "GlobalScalarData1";
    readVectorDataName      = "GlobalVectorData2";
    readScalarDataName      = "GlobalScalarData2";
    writeVectorValue        = 1;
    writeScalarValue        = 11;
    expectedReadVectorValue = 2;
    expectedReadScalarValue = 22;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    writeVectorDataName     = "GlobalVectorData2";
    writeScalarDataName     = "GlobalScalarData2";
    readVectorDataName      = "GlobalVectorData1";
    readScalarDataName      = "GlobalScalarData1";
    writeVectorValue        = 2;
    writeScalarValue        = 22;
    expectedReadVectorValue = 1;
    expectedReadScalarValue = 11;
  }
  int writeVectorDataID = couplingInterface.getGlobalDataID(writeVectorDataName);
  int readVectorDataID  = couplingInterface.getGlobalDataID(readVectorDataName);
  int writeScalarDataID = couplingInterface.getGlobalDataID(writeScalarDataName);
  int readScalarDataID  = couplingInterface.getGlobalDataID(readScalarDataName);

  double              dt = 0;
  std::vector<double> writeVectorData(dimensions, writeVectorValue);
  double              writeScalarData = writeScalarValue;
  std::vector<double> readVectorData(dimensions, -1);
  double              readScalarData = -1;

  dt = couplingInterface.initialize();

  while (couplingInterface.isCouplingOngoing()) {
    if (couplingInterface.requiresWritingCheckpoint()) {
    }
    // Write: from local data structure --> to precice buffer
    couplingInterface.writeGlobalVectorData(writeVectorDataID, writeVectorData.data());
    couplingInterface.writeGlobalScalarData(writeScalarDataID, writeScalarData);
    // Advance (exchange coupling data)
    dt = couplingInterface.advance(dt);
    // Read: from precice buffer --> to local data structure
    couplingInterface.readGlobalVectorData(readVectorDataID, readVectorData.data());
    couplingInterface.readGlobalScalarData(readScalarDataID, readScalarData);
    // Check read data
    BOOST_TEST(expectedReadVectorValue == readVectorData.at(0));
    BOOST_TEST(expectedReadVectorValue == readVectorData.at(1));
    BOOST_TEST(expectedReadScalarValue == readScalarData);
    if (couplingInterface.requiresReadingCheckpoint()) {
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // GlobalData

#endif // PRECICE_NO_MPI
