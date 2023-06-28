#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(GlobalData)

/**
 * @brief Test simple explicit coupled simulation with
 * global data initialization and without acceleration.
 * The second solver initializes the data of the first.
 *
 */
BOOST_AUTO_TEST_CASE(ExplicitInitializeData)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector3d;

  Participant participant(context.name, context.config(), 0, 1);

  if (context.isNamed("SolverOne")) {
    auto     dataAName = "GlobalData1";
    auto     dataBName = "GlobalData2";
    Vector3d valueDataB(0.0, 0.0, 0.0);
    Vector3d expectedDataB(2.0, 2.0, 2.0);

    participant.initialize();
    double maxDt = participant.getMaxTimeStepSize();
    participant.readGlobalData(dataBName, maxDt, valueDataB);
    BOOST_TEST(precice::testing::equals(expectedDataB, valueDataB));

    while (participant.isCouplingOngoing()) {
      Vector3d valueDataA(1.0, 1.0, 1.0);
      participant.writeGlobalData(dataAName, valueDataA);
      participant.advance(maxDt);
      maxDt = participant.getMaxTimeStepSize();
      participant.readGlobalData(dataBName, maxDt, valueDataB);
      Vector3d expectedDataB(2.5, 2.5, 2.5);
      BOOST_TEST(precice::testing::equals(expectedDataB, valueDataB));
    }
    participant.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    BOOST_REQUIRE(participant.requiresInitialData());
    auto dataAName = "GlobalData1";
    auto dataBName = "GlobalData2";

    Vector3d valueDataB(2.0, 2.0, 2.0);
    participant.writeGlobalData(dataBName, valueDataB);
    //tell preCICE that data has been written and call initializeData
    participant.initialize();
    double maxDt = participant.getMaxTimeStepSize();

    Vector3d valueDataA;
    participant.readGlobalData(dataAName, maxDt, valueDataA);
    Vector3d expectedDataA(1.0, 1.0, 1.0);
    BOOST_TEST(precice::testing::equals(valueDataA, expectedDataA));

    while (participant.isCouplingOngoing()) {
      Vector3d valueDataB(2.5, 2.5, 2.5);
      participant.writeGlobalData(dataBName, valueDataB);
      participant.advance(maxDt);
      maxDt = participant.getMaxTimeStepSize();
      participant.readGlobalData(dataAName, maxDt, valueDataA);
      BOOST_TEST(valueDataA == expectedDataA);
    }
    participant.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // GlobalData

#endif // PRECICE_NO_MPI
