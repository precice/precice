#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(InitializeData)

/**
 * @brief The second solver initializes the data of the first.
 *
 * A mapping is employed for the second solver, i.e., at the end of
 * initialize(), the mapping needs to be invoked.
 */
BOOST_AUTO_TEST_CASE(Explicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto   meshName   = "MeshOne";
    double v[]        = {1.0, 2.0, 3.0};
    auto   vid        = cplInterface.setMeshVertex(meshName, v);
    auto   dataAID    = "DataOne";
    auto   dataBID    = "DataTwo";
    double valueDataB = 0.0;
    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    cplInterface.readData(meshName, dataBID, {&vid, 1}, maxDt, {&valueDataB, 1});
    BOOST_TEST(2.0 == valueDataB);
    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(1.0, 1.0, 1.0);
      cplInterface.writeData(meshName, dataAID, {&vid, 1}, valueDataA);
      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();
      cplInterface.readData(meshName, dataBID, {&vid, 1}, maxDt, {&valueDataB, 1});
      BOOST_TEST(2.5 == valueDataB);
    }
    cplInterface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto     meshName = "MeshTwo";
    Vector3d pos      = Vector3d::Zero();
    auto     vid      = cplInterface.setMeshVertex(meshName, pos);

    BOOST_REQUIRE(cplInterface.requiresInitialData());
    auto   dataAID      = "DataOne";
    auto   dataBID      = "DataTwo";
    double valueDataB[] = {2.0};
    cplInterface.writeData(meshName, dataBID, {&vid, 1}, valueDataB);
    //tell preCICE that data has been written and call initializeData
    cplInterface.initialize();
    double   maxDt = cplInterface.getMaxTimeStepSize();
    Vector3d valueDataA;
    cplInterface.readData(meshName, dataAID, {&vid, 1}, maxDt, valueDataA);
    Vector3d expected(1.0, 1.0, 1.0);
    BOOST_TEST(valueDataA == expected);
    while (cplInterface.isCouplingOngoing()) {
      double valueDataB[] = {2.5};
      cplInterface.writeData(meshName, dataBID, {&vid, 1}, valueDataB);
      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();
      cplInterface.readData(meshName, dataAID, {&vid, 1}, maxDt, valueDataA);
      BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // InitializeData

#endif // PRECICE_NO_MPI
