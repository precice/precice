#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>
#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(GeomultiscaleMapping)
BOOST_AUTO_TEST_CASE(DimensionsConfiguration)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

  /*  In this test case, SolverOne is the 1D code (despite having 3D vertices, due to current shortcomings)
      and we're testing the AxialGeoMultiscaleMapping feature with a parabolic inlet profile at the (downstream) 3D inlet.
  */

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto     meshName = "MeshOne";
    Vector3d posOne   = Vector3d::Constant(0.0);
    auto     vid      = cplInterface.setMeshVertex(meshName, posOne);
    auto     dataAID  = "DataOne";
    auto     dataBID  = "DataTwo";

    Vector3d valueDataB;

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    cplInterface.readData(meshName, dataBID, {&vid, 1}, maxDt, valueDataB);
    Vector3d expected(0.0, 0.0, 8.0);
    // BOOST_TEST(valueDataB == expected);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(0.0, 0.0, 8.0);
      cplInterface.writeData(meshName, dataAID, {&vid, 1}, valueDataA);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName, dataBID, {&vid, 1}, maxDt, valueDataB);
      expected << 0.0, 0.0, 8.0;
      // BOOST_TEST(valueDataB == expected);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto     meshName = "MeshTwo";
    Vector3d pos      = Vector3d::Constant(0.0);
    auto     vid      = cplInterface.setMeshVertex(meshName, pos);

    auto dataAID = "DataOne";
    auto dataBID = "DataTwo";

    BOOST_REQUIRE(cplInterface.requiresInitialData());
    Vector3d valueDataB(0.0, 0.0, 4.0);
    cplInterface.writeData(meshName, dataBID, {&vid, 1}, valueDataB);

    //tell preCICE that data has been written and call initialize
    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    Vector3d valueDataA;
    cplInterface.readData(meshName, dataAID, {&vid, 1}, maxDt, valueDataA);
    Vector3d expected(0.0, 0.0, 16.0);
    // BOOST_TEST(valueDataA == expected);

    while (cplInterface.isCouplingOngoing()) {

      valueDataB << 0.0, 0.0, 4.0;
      cplInterface.writeData(meshName, dataBID, {&vid, 1}, valueDataB);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();
      cplInterface.readData(meshName, dataAID, {&vid, 1}, maxDt, valueDataA);
      // BOOST_TEST(valueDataA == expected);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // GeomultiscaleMapping

#endif // PRECICE_NO_MPI
