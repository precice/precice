#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>
#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(GeometricMultiscale)
BOOST_AUTO_TEST_CASE(AxialGeoMultiscale)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))

  /*  In this test case, SolverOne is the 1D code (despite having 3D vertices, due to current shortcomings)
      and we're testing the AxialGeoMultiscaleMapping feature with a parabolic inlet profile at the (downstream) 3D inlet.

      We are assuming that the exchanged data follows a parabolic velocity profile, which means that a 3D
      maximum velocity value is mapped to a 1D average velocity value. Therefore, the expected value in the
      test is different than the one set on the other side.
  */

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), 0, 1);
  if (context.isNamed("SolverOne")) {
    auto     meshName  = "MeshOne";
    Vector3d posOne    = Vector3d::Constant(0.0);
    auto     vid       = cplInterface.setMeshVertex(meshName, posOne);
    auto     dataAName = "DataOne";
    auto     dataBName = "DataTwo";

    Vector3d valueDataB;

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();
    cplInterface.readData(meshName, dataBName, {&vid, 1}, maxDt, valueDataB);
    Vector3d expectedDataB(0.0, 0.0, 8.0); // Read, averaged (collect) from one vertex with value 8.0
    BOOST_TEST(valueDataB == expectedDataB);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(0.0, 0.0, 8.0); // Written, to be scaled up (spread) to 16.0, following a parabolic profile
      cplInterface.writeData(meshName, dataAName, {&vid, 1}, valueDataA);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName, dataBName, {&vid, 1}, maxDt, valueDataB);
      BOOST_TEST(valueDataB == expectedDataB);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    auto meshName = "MeshTwo";
    // This test is only defining one vertex in the 3D mesh, and the collect operation
    // is averaging. Therefore, the mapped value will always be the same, since we
    // are only defining one point.
    Vector3d pos       = Vector3d::Constant(0.0);
    auto     vid       = cplInterface.setMeshVertex(meshName, pos);
    auto     dataAName = "DataOne";
    auto     dataBName = "DataTwo";

    BOOST_REQUIRE(cplInterface.requiresInitialData());
    Vector3d valueDataB(0.0, 0.0, 8.0); // Written, to be averaged (collect) to 8.0
    cplInterface.writeData(meshName, dataBName, {&vid, 1}, valueDataB);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    Vector3d valueDataA;
    cplInterface.readData(meshName, dataAName, {&vid, 1}, maxDt, valueDataA);
    Vector3d expectedDataA(0.0, 0.0, 16.0); // Read, scaled up (spread) from 8.0
    BOOST_TEST(valueDataA == expectedDataA);

    while (cplInterface.isCouplingOngoing()) {

      valueDataB << 0.0, 0.0, 8.0; // Written, to be averaged (collect) to 8.0
      cplInterface.writeData(meshName, dataBName, {&vid, 1}, valueDataB);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();
      cplInterface.readData(meshName, dataAName, {&vid, 1}, maxDt, valueDataA);
      BOOST_TEST(valueDataA == expectedDataA);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // GeomultiscaleMapping

#endif // PRECICE_NO_MPI
