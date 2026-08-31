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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleScalarUniform1D3D_Square)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("Fluid1D")) {
    const auto meshName1D   = "Mesh1D";
    const auto dataWriteSca = "PressureLike";
    const auto dataReadVec  = "VelocityLike";

    Vector3d pos(0.0, 0.0, 0.0);
    auto     vid = cplInterface.setMeshVertex(meshName1D, pos);

    Vector3d valueReadVec = Vector3d::Zero();

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    cplInterface.readData(meshName1D, dataReadVec, {&vid, 1}, maxDt, valueReadVec);
    Vector3d expectedVec(0.0, 0.0, 6.0);
    BOOST_TEST(valueReadVec == expectedVec);

    while (cplInterface.isCouplingOngoing()) {
      double valueWrite = 8.0;
      cplInterface.writeData(meshName1D, dataWriteSca, {&vid, 1}, {&valueWrite, 1});

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName1D, dataReadVec, {&vid, 1}, maxDt, valueReadVec);
      BOOST_TEST(valueReadVec == expectedVec);
    }

    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("Fluid3D"));

    const auto meshName3D   = "Mesh3D";
    const auto dataReadSca  = "PressureLike";
    const auto dataWriteVec = "VelocityLike";

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

    Eigen::VectorXd vInit(12);
    vInit << 0.0, 0.0, 8.0,
        0.0, 0.0, 4.0,
        0.0, 0.0, 6.0,
        0.0, 0.0, 6.0; // average = 6

    cplInterface.writeData(meshName3D, dataWriteVec, vids3D, vInit);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    double pRead0 = 0.0, pRead1 = 0.0, pRead2 = 0.0, pRead3 = 0.0;

    cplInterface.readData(meshName3D, dataReadSca, {&vid0, 1}, maxDt, {&pRead0, 1});
    cplInterface.readData(meshName3D, dataReadSca, {&vid1, 1}, maxDt, {&pRead1, 1});
    cplInterface.readData(meshName3D, dataReadSca, {&vid2, 1}, maxDt, {&pRead2, 1});
    cplInterface.readData(meshName3D, dataReadSca, {&vid3, 1}, maxDt, {&pRead3, 1});

    BOOST_TEST(pRead0 == 8.0);
    BOOST_TEST(pRead1 == 8.0);
    BOOST_TEST(pRead2 == 8.0);
    BOOST_TEST(pRead3 == 8.0);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName3D, dataWriteVec, vids3D, vInit);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName3D, dataReadSca, {&vid0, 1}, maxDt, {&pRead0, 1});
      cplInterface.readData(meshName3D, dataReadSca, {&vid1, 1}, maxDt, {&pRead1, 1});
      cplInterface.readData(meshName3D, dataReadSca, {&vid2, 1}, maxDt, {&pRead2, 1});
      cplInterface.readData(meshName3D, dataReadSca, {&vid3, 1}, maxDt, {&pRead3, 1});

      BOOST_TEST(pRead0 == 8.0);
      BOOST_TEST(pRead1 == 8.0);
      BOOST_TEST(pRead2 == 8.0);
      BOOST_TEST(pRead3 == 8.0);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
