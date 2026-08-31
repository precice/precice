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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleVectorUniform1D3D_Square)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("Fluid1D")) {
    const auto meshName1D   = "Mesh1D";
    const auto dataWriteVec = "VelocityLike";
    const auto dataReadSca  = "PressureLike";

    Vector3d pos(0.0, 0.0, 0.0);
    auto     vid = cplInterface.setMeshVertex(meshName1D, pos);

    double valueReadSca = 0.0;

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    cplInterface.readData(meshName1D, dataReadSca, {&vid, 1}, maxDt, {&valueReadSca, 1});
    const double expectedSca = 6.0;
    BOOST_TEST(valueReadSca == expectedSca);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueWrite(0.0, 0.0, 8.0);
      cplInterface.writeData(meshName1D, dataWriteVec, {&vid, 1}, valueWrite);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName1D, dataReadSca, {&vid, 1}, maxDt, {&valueReadSca, 1});
      BOOST_TEST(valueReadSca == expectedSca);
    }

    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("Fluid3D"));

    const auto meshName3D   = "Mesh3D";
    const auto dataReadVec  = "VelocityLike";
    const auto dataWriteSca = "PressureLike";

    Vector3d p0(0.0, 0.0, 0.0);
    Vector3d p1(1.0, 0.0, 0.0);
    Vector3d p2(0.0, 1.0, 0.0);
    Vector3d p3(1.0, 1.0, 0.0);

    auto vid0 = cplInterface.setMeshVertex(meshName3D, p0);
    auto vid1 = cplInterface.setMeshVertex(meshName3D, p1);
    auto vid2 = cplInterface.setMeshVertex(meshName3D, p2);
    auto vid3 = cplInterface.setMeshVertex(meshName3D, p3);

    std::array<int, 4> vids3D = {vid0, vid1, vid2, vid3};

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    Eigen::VectorXd pInit(4);
    pInit << 8.0, 4.0, 6.0, 6.0; // average = 6
    cplInterface.writeData(meshName3D, dataWriteSca, vids3D, pInit);

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

    Vector3d exp(0.0, 0.0, 8.0);

    BOOST_TEST(vOut0 == exp);
    BOOST_TEST(vOut1 == exp);
    BOOST_TEST(vOut2 == exp);
    BOOST_TEST(vOut3 == exp);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName3D, dataWriteSca, vids3D, pInit);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName3D, dataReadVec, {&vid0, 1}, maxDt, vOut0);
      cplInterface.readData(meshName3D, dataReadVec, {&vid1, 1}, maxDt, vOut1);
      cplInterface.readData(meshName3D, dataReadVec, {&vid2, 1}, maxDt, vOut2);
      cplInterface.readData(meshName3D, dataReadVec, {&vid3, 1}, maxDt, vOut3);

      BOOST_TEST(vOut0 == exp);
      BOOST_TEST(vOut1 == exp);
      BOOST_TEST(vOut2 == exp);
      BOOST_TEST(vOut3 == exp);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
