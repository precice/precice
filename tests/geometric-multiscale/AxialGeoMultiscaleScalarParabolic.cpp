#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>
#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(GeometricMultiscale)
PRECICE_TEST_SETUP("Fluid1D"_on(1_rank), "Fluid3D"_on(1_rank))
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleScalarParabolic)
{
  PRECICE_TEST();

  /*  Reverse case:
      - Fluid1D writes PressureLike (to be SPREAD with in a parabolic way across the vertices)
      - Fluid1D reads VelocityLike (COLLECT/averaged from 3D)
      This mirrors the non-reverse test by swapping which field is spread/collected.
  */

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);
  if (context.isNamed("Fluid1D")) {
    auto     meshName  = "Mesh1D";
    Vector3d posOne    = Vector3d::Constant(0.0);
    auto     vid       = cplInterface.setMeshVertex(meshName, posOne);
    auto     dataAName = "PressureLike";
    auto     dataBName = "VelocityLike";

    Vector3d valueDataB;
    double   valueDataA;

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // First read: VelocityLike collected from 3D. With one 3D vertex at (0,0,8), average stays 8.
    cplInterface.readData(meshName, dataBName, {&vid, 1}, maxDt, valueDataB);
    Vector3d expectedDataB(0.0, 0.0, 8.0);
    BOOST_TEST(valueDataB == expectedDataB);

    while (cplInterface.isCouplingOngoing()) {
      // Write scalar PressureLike = 8.0 on 1D -> SPREAD to 3D (parabolic, scaled up to 16.0 at 3D centerline)
      valueDataA = 8.0;
      cplInterface.writeData(meshName, dataAName, {&vid, 1}, {&valueDataA, 1});

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      // Read VelocityLike again (collected from 3D); still (0,0,8) with a single 3D vertex providing 8
      cplInterface.readData(meshName, dataBName, {&vid, 1}, maxDt, valueDataB);
      BOOST_TEST(valueDataB == expectedDataB);
    }
    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("Fluid3D"));
    auto meshName = "Mesh3D";
    // This test is only defining one vertex in the 3D mesh, and the collect operation
    // is averaging. Therefore, the mapped value will always be the same, since we
    // are only defining one point.
    Vector3d pos       = Vector3d::Constant(0.0);
    auto     vid       = cplInterface.setMeshVertex(meshName, pos);
    auto     dataAName = "PressureLike";
    auto     dataBName = "VelocityLike";

    BOOST_REQUIRE(cplInterface.requiresInitialData());
    Vector3d valueDataB(0.0, 0.0, 8.0); // Written, to be averaged (collect) to 8.0
    cplInterface.writeData(meshName, dataBName, {&vid, 1}, valueDataB);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    double valueDataA;
    cplInterface.readData(meshName, dataAName, {&vid, 1}, maxDt, {&valueDataA, 1});
    double expectedDataA = 16.0;
    BOOST_TEST(valueDataA == expectedDataA);

    while (cplInterface.isCouplingOngoing()) {

      valueDataB << 0.0, 0.0, 8.0; // Written, to be averaged (collect) to 8.0
      cplInterface.writeData(meshName, dataBName, {&vid, 1}, valueDataB);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();
      cplInterface.readData(meshName, dataAName, {&vid, 1}, maxDt, {&valueDataA, 1});
      BOOST_TEST(valueDataA == expectedDataA);
    }
    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // GeomultiscaleMapping

#endif // PRECICE_NO_MPI
