#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(GlobalData)
// Test case to send global (meshless) data from one participant to another.
// The test case here is the most basic variant in order
// use such a feature. SolverOne defines the data and SolverTwo receives it.
BOOST_AUTO_TEST_CASE(Explicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  precice::SolverInterface couplingInterface(context.name, context.config(), 0, 1);
  const int                dimensions = 2;

  if (context.isNamed("SolverOne")) {
    const std::string globalDataName       = "GlobalData1";
    const std::string globalVectorDataName = "GlobalVectorData";
    couplingInterface.initialize();
    double dt = couplingInterface.getMaxTimeStepSize();
    // Some dummy writeData
    double              writeGlobalData{5};
    std::vector<double> writeGlobalVectorData(dimensions, 50.5);

    while (couplingInterface.isCouplingOngoing()) {
      // Write data to be sent to SolverTwo to buffer
      couplingInterface.writeGlobalScalarData(globalDataName, writeGlobalData);
      couplingInterface.writeGlobalVectorData(globalVectorDataName, writeGlobalVectorData.data());
      // send data
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
      // change reference data for next check
      writeGlobalData++;
      for (auto &elem : writeGlobalVectorData) {
        elem++;
      }
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    const std::string globalDataName       = "GlobalData1";
    const std::string globalVectorDataName = "GlobalVectorData";

    // Allocate data to read
    double              readGlobalData;
    std::vector<double> readGlobalVectorData(dimensions, -1);

    // Initialize
    double              expectedGlobalData{5};
    std::vector<double> expectedGlobalVectorData(dimensions, 50.5);
    couplingInterface.initialize(); // For serial-explicit, first communication happens here
    double dt = couplingInterface.getMaxTimeStepSize();

    while (couplingInterface.isCouplingOngoing()) {
      // read received data from buffer
      couplingInterface.readGlobalScalarData(globalDataName, readGlobalData);
      couplingInterface.readGlobalVectorData(globalVectorDataName, readGlobalVectorData.data());
      // check if received data is correct
      BOOST_TEST(precice::testing::equals(expectedGlobalData, readGlobalData));
      BOOST_TEST(precice::testing::equals(expectedGlobalVectorData, readGlobalVectorData));
      // receive next data
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
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
