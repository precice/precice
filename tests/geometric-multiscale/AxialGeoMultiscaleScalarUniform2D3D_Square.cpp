#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>
#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(GeometricMultiscale)

PRECICE_TEST_SETUP("Fluid2D"_on(1_rank), "Fluid3D"_on(1_rank))
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleScalarUniform2D3D_Square)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("Fluid2D")) {
    const auto meshName2D   = "Mesh2D";
    const auto dataWriteSca = "PressureLike";
    const auto dataReadVec  = "VelocityLike";

    Vector3d v0(0.0, -0.5, 0.0);
    Vector3d v1(0.0, 0.5, 0.0);

    auto vid0 = cplInterface.setMeshVertex(meshName2D, v0);
    auto vid1 = cplInterface.setMeshVertex(meshName2D, v1);

    std::array<int, 2> vids2D = {vid0, vid1};

    Vector3d valueReadVec = Vector3d::Zero();

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    cplInterface.readData(meshName2D, dataReadVec, {&vid0, 1}, maxDt, valueReadVec);
    Vector3d expectedVec(0.0, 0.0, 6.0);
    BOOST_TEST(valueReadVec == expectedVec);

    while (cplInterface.isCouplingOngoing()) {
      Eigen::VectorXd pIn(2);
      pIn << 10.0, 20.0;
      cplInterface.writeData(meshName2D, dataWriteSca, vids2D, pIn);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName2D, dataReadVec, {&vid0, 1}, maxDt, valueReadVec);
      BOOST_TEST(valueReadVec == expectedVec);
    }

    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("Fluid3D"));

    const auto meshName3D   = "Mesh3D";
    const auto dataReadSca  = "PressureLike";
    const auto dataWriteVec = "VelocityLike";

    Vector3d w0(0.0, -0.5, 0.0); // nearest v0
    Vector3d w1(0.5, -0.5, 0.0); // nearest v0
    Vector3d w2(0.0, 0.5, 0.0);  // nearest v1
    Vector3d w3(0.5, 0.5, 0.0);  // nearest v1

    auto wid0 = cplInterface.setMeshVertex(meshName3D, w0);
    auto wid1 = cplInterface.setMeshVertex(meshName3D, w1);
    auto wid2 = cplInterface.setMeshVertex(meshName3D, w2);
    auto wid3 = cplInterface.setMeshVertex(meshName3D, w3);

    std::array<int, 4> vids3D = {wid0, wid1, wid2, wid3};

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    Eigen::VectorXd vInit(12);
    vInit << 0.0, 0.0, 8.0,
        0.0, 0.0, 4.0,
        0.0, 0.0, 6.0,
        0.0, 0.0, 6.0;

    cplInterface.writeData(meshName3D, dataWriteVec, vids3D, vInit);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    double p0 = 0.0, p1 = 0.0, p2 = 0.0, p3 = 0.0;

    cplInterface.readData(meshName3D, dataReadSca, {&wid0, 1}, maxDt, {&p0, 1});
    cplInterface.readData(meshName3D, dataReadSca, {&wid1, 1}, maxDt, {&p1, 1});
    cplInterface.readData(meshName3D, dataReadSca, {&wid2, 1}, maxDt, {&p2, 1});
    cplInterface.readData(meshName3D, dataReadSca, {&wid3, 1}, maxDt, {&p3, 1});

    BOOST_TEST(p0 == 10.0);
    BOOST_TEST(p1 == 10.0);
    BOOST_TEST(p2 == 20.0);
    BOOST_TEST(p3 == 20.0);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName3D, dataWriteVec, vids3D, vInit);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName3D, dataReadSca, {&wid0, 1}, maxDt, {&p0, 1});
      cplInterface.readData(meshName3D, dataReadSca, {&wid1, 1}, maxDt, {&p1, 1});
      cplInterface.readData(meshName3D, dataReadSca, {&wid2, 1}, maxDt, {&p2, 1});
      cplInterface.readData(meshName3D, dataReadSca, {&wid3, 1}, maxDt, {&p3, 1});

      BOOST_TEST(p0 == 10.0);
      BOOST_TEST(p1 == 10.0);
      BOOST_TEST(p2 == 20.0);
      BOOST_TEST(p3 == 20.0);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
