#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>
#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(GeometricMultiscale)

PRECICE_TEST_SETUP("Fluid1D"_on(1_rank), "Fluid2D"_on(1_rank))
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleScalarParabolic1D2D)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("Fluid1D")) {
    const auto meshName = "Mesh1D";
    Vector3d   posOne   = Vector3d::Constant(0.0);
    auto       vid      = cplInterface.setMeshVertex(meshName, posOne);

    const auto dataAName = "PressureLike";
    const auto dataBName = "VelocityLike";

    Vector3d valueDataB = Vector3d::Zero();
    double   valueDataA = 0.0;

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    cplInterface.readData(meshName, dataBName, {&vid, 1}, maxDt, valueDataB);
    Vector3d expectedDataB(0.0, 0.0, 6.0);
    BOOST_TEST(valueDataB == expectedDataB);

    while (cplInterface.isCouplingOngoing()) {
      valueDataA = 8.0;
      cplInterface.writeData(meshName, dataAName, {&vid, 1}, {&valueDataA, 1});

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName, dataBName, {&vid, 1}, maxDt, valueDataB);
      BOOST_TEST(valueDataB == expectedDataB);
    }

    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("Fluid2D"));

    const auto meshName  = "Mesh2D";
    const auto dataAName = "PressureLike";
    const auto dataBName = "VelocityLike";

    Vector3d p0(0.0, 0.0, 0.0);
    Vector3d p1(0.0, 0.5, 0.0);

    auto vid0 = cplInterface.setMeshVertex(meshName, p0);
    auto vid1 = cplInterface.setMeshVertex(meshName, p1);

    std::array<int, 2> vids = {vid0, vid1};

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    Vector3d v0(0.0, 0.0, 8.0);
    Vector3d v1(0.0, 0.0, 4.0);

    Eigen::VectorXd initValues(6);
    initValues << v0(0), v0(1), v0(2),
        v1(0), v1(1), v1(2);

    cplInterface.writeData(meshName, dataBName, vids, initValues);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeData(meshName, dataBName, vids, initValues);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      double pOut0 = 0.0;
      double pOut1 = 0.0;

      cplInterface.readData(meshName, dataAName, {&vid0, 1}, maxDt, {&pOut0, 1});
      cplInterface.readData(meshName, dataAName, {&vid1, 1}, maxDt, {&pOut1, 1});

      constexpr double factor    = 1.5;
      const double     expected0 = factor * 8.0 * (1.0 - 0.0 * 0.0);
      const double     expected1 = factor * 8.0 * (1.0 - 0.5 * 0.5);

      BOOST_TEST(pOut0 == expected0);
      BOOST_TEST(pOut1 == expected1);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
