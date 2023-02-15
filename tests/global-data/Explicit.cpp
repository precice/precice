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

  // Set up Solverinterface
  precice::SolverInterface couplingInterface(context.name, context.config(), 0, 1);
  BOOST_TEST(couplingInterface.getDimensions() == 2);

  if (context.isNamed("SolverOne")) {
    const int globalDataID = couplingInterface.getGlobalDataID("GlobalData1");
    double    dt           = couplingInterface.initialize();
    // Some dummy writeData
    double writeGlobalData{5};

    while (couplingInterface.isCouplingOngoing()) {
      // Write data
      couplingInterface.writeGlobalScalarData(globalDataID, writeGlobalData);
      dt = couplingInterface.advance(dt);
    }

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    // Query IDs
    const int globalDataID = couplingInterface.getGlobalDataID("GlobalData1");

    // Allocate data to read
    double readGlobalData;

    // Initialize
    double dt = couplingInterface.initialize();
    while (couplingInterface.isCouplingOngoing()) {
      dt = couplingInterface.advance(dt);
      couplingInterface.readGlobalScalarData(globalDataID, readGlobalData);
      // Expected data according to the writeData
      double expectedGlobalData{5};
      BOOST_TEST(precice::testing::equals(expectedGlobalData, readGlobalData));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // GlobalData

#endif // PRECICE_NO_MPI
