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
    const std::string globalScalarDataName = "GlobalData1";
    const std::string globalVectorDataName = "GlobalVectorData";
    couplingInterface.initialize();
    double dt = couplingInterface.getMaxTimeStepSize();
    // Some dummy writeData
    double              writeGlobalScalarData{5};
    std::vector<double> writeGlobalVectorData(dimensions, 50.5);

    while (couplingInterface.isCouplingOngoing()) {
      // Write data to be sent to SolverTwo to buffer
      couplingInterface.writeGlobalData(globalScalarDataName, {&writeGlobalScalarData, 1});
      couplingInterface.writeGlobalData(globalVectorDataName, writeGlobalVectorData);
      // send data
      couplingInterface.advance(dt);
      dt = couplingInterface.getMaxTimeStepSize();
      // change reference data for next check
      writeGlobalScalarData++;
      for (auto &elem : writeGlobalVectorData) {
        elem++;
      }
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    const std::string globalScalarDataName = "GlobalData1";
    const std::string globalVectorDataName = "GlobalVectorData";

    // Allocate data to read
    double              readGlobalScalarData;
    std::vector<double> readGlobalVectorData(dimensions, -1);

    // Initialize
    double              expectedGlobalData{5};
    std::vector<double> expectedGlobalVectorData(dimensions, 50.5);
    couplingInterface.initialize(); // For serial-explicit, first communication happens here
    double dt = couplingInterface.getMaxTimeStepSize();

    while (couplingInterface.isCouplingOngoing()) {
      // read received data from buffer
      couplingInterface.readGlobalScalarData(globalScalarDataName, readGlobalScalarData);
      couplingInterface.readGlobalVectorData(globalVectorDataName, readGlobalVectorData.data());
      // check if received data is correct
      BOOST_TEST(precice::testing::equals(expectedGlobalData, readGlobalScalarData));
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
