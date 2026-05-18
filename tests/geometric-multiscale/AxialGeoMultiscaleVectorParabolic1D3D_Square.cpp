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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleVectorParabolic1D3D_Square)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("Fluid1D")) {
    const auto meshNameVec = "Mesh1D";
    Vector3d   posOne      = Vector3d::Constant(0.0);
    auto       vid         = cplInterface.setMeshVertex(meshNameVec, posOne);

    const auto dataWriteVec = "VelocityLike";
    const auto dataReadSca  = "PressureLike";

    double valueReadScalar = 0.0;

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    cplInterface.readData(meshNameVec, dataReadSca, {&vid, 1}, maxDt, {&valueReadScalar, 1});
    const double expectedScalar = 6.0;
    BOOST_TEST(valueReadScalar == expectedScalar);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueWriteVec(0.0, 0.0, 8.0);
      cplInterface.writeData(meshNameVec, dataWriteVec, {&vid, 1}, valueWriteVec);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshNameVec, dataReadSca, {&vid, 1}, maxDt, {&valueReadScalar, 1});
      BOOST_TEST(valueReadScalar == expectedScalar);
    }

    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("Fluid3D"));

    const auto meshName3D   = "Mesh3D";
    const auto dataReadVec  = "VelocityLike";
    const auto dataWriteSca = "PressureLike";

    Vector3d p0(0.0, 0.0, 0.0); // center
    Vector3d p1(1.0, 0.0, 0.0); // edge (s1=1,s2=0) -> 0
    Vector3d p2(0.0, 0.5, 0.0); // mid  (s1=0,s2=0.5)
    Vector3d p3(0.5, 0.5, 0.0); // corner (s1=0.5,s2=0.5)

    auto vid0 = cplInterface.setMeshVertex(meshName3D, p0);
    auto vid1 = cplInterface.setMeshVertex(meshName3D, p1);
    auto vid2 = cplInterface.setMeshVertex(meshName3D, p2);
    auto vid3 = cplInterface.setMeshVertex(meshName3D, p3);

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    std::array<int, 4> vids = {vid0, vid1, vid2, vid3};

    Eigen::VectorXd initPressure(4);
    initPressure << 8.0, 4.0, 6.0, 6.0; // average = 6.0
    cplInterface.writeData(meshName3D, dataWriteSca, vids, initPressure);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    Vector3d vOut0 = Vector3d::Zero();
    Vector3d vOut1 = Vector3d::Zero();
    Vector3d vOut2 = Vector3d::Zero();
    Vector3d vOut3 = Vector3d::Zero();

    cplInterface.readData(meshName3D, dataReadVec, {&vid0, 1}, maxDt, vOut0);
    cplInterface.readData(meshName3D, dataReadVec, {&vid1, 1}, maxDt, vOut1);
    cplInterface.readData(meshName3D, dataReadVec, {&vid2, 1}, maxDt, vOut2);
    cplInterface.readData(meshName3D, dataReadVec, {&vid3, 1}, maxDt, vOut3);

    constexpr double factor = 2.096;
    constexpr double m      = 0.879;

    const double u = 8.0;

    const double b00 = 1.0;
    const double b05 = 1.0 - 0.5 * 0.5; // 0.75

    const double expected0 = factor * u * std::pow(b00, m) * std::pow(b00, m);
    const double expected1 = 0.0;
    const double expected2 = factor * u * std::pow(b00, m) * std::pow(b05, m);
    const double expected3 = factor * u * std::pow(b05, m) * std::pow(b05, m);

    Vector3d exp0(0.0, 0.0, expected0);
    Vector3d exp1(0.0, 0.0, expected1);
    Vector3d exp2(0.0, 0.0, expected2);
    Vector3d exp3(0.0, 0.0, expected3);

    BOOST_TEST(vOut0 == exp0);
    BOOST_TEST(vOut1 == exp1);
    BOOST_TEST(vOut2 == exp2);
    BOOST_TEST(vOut3 == exp3);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName3D, dataWriteSca, vids, initPressure);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName3D, dataReadVec, {&vid0, 1}, maxDt, vOut0);
      cplInterface.readData(meshName3D, dataReadVec, {&vid1, 1}, maxDt, vOut1);
      cplInterface.readData(meshName3D, dataReadVec, {&vid2, 1}, maxDt, vOut2);
      cplInterface.readData(meshName3D, dataReadVec, {&vid3, 1}, maxDt, vOut3);

      BOOST_TEST(vOut0 == exp0);
      BOOST_TEST(vOut1 == exp1);
      BOOST_TEST(vOut2 == exp2);
      BOOST_TEST(vOut3 == exp3);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
