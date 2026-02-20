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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleScalarParabolic1D3D_Square)
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
    Vector3d expectedDataB(0.0, 0.0, 8.0);
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
    BOOST_TEST(context.isNamed("Fluid3D"));

    const auto meshName  = "Mesh3D";
    const auto dataAName = "PressureLike";
    const auto dataBName = "VelocityLike";

    Vector3d p0(0.0, 0.0, 0.0);
    Vector3d p1(1.0, 0.0, 0.0);
    Vector3d p2(0.0, 0.5, 0.0);
    Vector3d p3(0.5, 0.5, 0.0);

    auto vid0 = cplInterface.setMeshVertex(meshName, p0);
    auto vid1 = cplInterface.setMeshVertex(meshName, p1);
    auto vid2 = cplInterface.setMeshVertex(meshName, p2);
    auto vid3 = cplInterface.setMeshVertex(meshName, p3);

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    Vector3d valueDataB(0.0, 0.0, 8.0);
    cplInterface.writeData(meshName, dataBName, {&vid0, 1}, valueDataB);
    cplInterface.writeData(meshName, dataBName, {&vid1, 1}, valueDataB);
    cplInterface.writeData(meshName, dataBName, {&vid2, 1}, valueDataB);
    cplInterface.writeData(meshName, dataBName, {&vid3, 1}, valueDataB);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    double pOut0 = 0.0, pOut1 = 0.0, pOut2 = 0.0, pOut3 = 0.0;
    cplInterface.readData(meshName, dataAName, {&vid0, 1}, maxDt, {&pOut0, 1});
    cplInterface.readData(meshName, dataAName, {&vid1, 1}, maxDt, {&pOut1, 1});
    cplInterface.readData(meshName, dataAName, {&vid2, 1}, maxDt, {&pOut2, 1});
    cplInterface.readData(meshName, dataAName, {&vid3, 1}, maxDt, {&pOut3, 1});

    constexpr double factor = 2.096;
    constexpr double m      = 0.879;

    const double b00 = 1.0;
    const double b10 = 0.0;
    const double b05 = 1.0 - 0.5 * 0.5;

    const double expected0 = factor * 8.0 * std::pow(b00, m) * std::pow(b00, m);
    const double expected1 = 0.0;
    const double expected2 = factor * 8.0 * std::pow(b00, m) * std::pow(b05, m);
    const double expected3 = factor * 8.0 * std::pow(b05, m) * std::pow(b05, m);

    BOOST_TEST(pOut0 == expected0);
    BOOST_TEST(pOut1 == expected1);
    BOOST_TEST(pOut2 == expected2);
    BOOST_TEST(pOut3 == expected3);

    while (cplInterface.isCouplingOngoing()) {

      valueDataB << 0.0, 0.0, 8.0;
      cplInterface.writeData(meshName, dataBName, {&vid0, 1}, valueDataB);
      cplInterface.writeData(meshName, dataBName, {&vid1, 1}, valueDataB);
      cplInterface.writeData(meshName, dataBName, {&vid2, 1}, valueDataB);
      cplInterface.writeData(meshName, dataBName, {&vid3, 1}, valueDataB);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName, dataAName, {&vid0, 1}, maxDt, {&pOut0, 1});
      cplInterface.readData(meshName, dataAName, {&vid1, 1}, maxDt, {&pOut1, 1});
      cplInterface.readData(meshName, dataAName, {&vid2, 1}, maxDt, {&pOut2, 1});
      cplInterface.readData(meshName, dataAName, {&vid3, 1}, maxDt, {&pOut3, 1});

      BOOST_TEST(pOut0 == expected0);
      BOOST_TEST(pOut1 == expected1);
      BOOST_TEST(pOut2 == expected2);
      BOOST_TEST(pOut3 == expected3);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
