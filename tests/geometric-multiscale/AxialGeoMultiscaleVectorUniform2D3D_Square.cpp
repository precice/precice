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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleVectorUniform2D3D_Square)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("Fluid2D")) {
    const auto meshName2D   = "Mesh2D";
    const auto dataWriteVec = "VelocityLike";
    const auto dataReadSca  = "PressureLike";

    Vector3d v0(0.0, -0.5, 0.0);
    Vector3d v1(0.0, 0.5, 0.0);

    auto vid0 = cplInterface.setMeshVertex(meshName2D, v0);
    auto vid1 = cplInterface.setMeshVertex(meshName2D, v1);

    std::array<int, 2> vids2D = {vid0, vid1};

    double valueReadSca = 0.0;

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    cplInterface.readData(meshName2D, dataReadSca, {&vid0, 1}, maxDt, {&valueReadSca, 1});
    const double expectedSca = 6.0;
    BOOST_TEST(valueReadSca == expectedSca);

    while (cplInterface.isCouplingOngoing()) {
      Eigen::VectorXd vIn(6);
      vIn << 0.0, 0.0, 10.0,
          0.0, 0.0, 20.0;
      cplInterface.writeData(meshName2D, dataWriteVec, {vids2D.data(), 2}, vIn);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName2D, dataReadSca, {&vid0, 1}, maxDt, {&valueReadSca, 1});
      BOOST_TEST(valueReadSca == expectedSca);
    }

    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("Fluid3D"));

    const auto meshName3D   = "Mesh3D";
    const auto dataReadVec  = "VelocityLike";
    const auto dataWriteSca = "PressureLike";

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

    Eigen::VectorXd pInit(4);
    pInit << 8.0, 4.0, 6.0, 6.0; // average = 6.0
    cplInterface.writeData(meshName3D, dataWriteSca, {vids3D.data(), 4}, pInit);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    Vector3d vOut0 = Vector3d::Zero();
    Vector3d vOut1 = Vector3d::Zero();
    Vector3d vOut2 = Vector3d::Zero();
    Vector3d vOut3 = Vector3d::Zero();

    cplInterface.readData(meshName3D, dataReadVec, {&wid0, 1}, maxDt, vOut0);
    cplInterface.readData(meshName3D, dataReadVec, {&wid1, 1}, maxDt, vOut1);
    cplInterface.readData(meshName3D, dataReadVec, {&wid2, 1}, maxDt, vOut2);
    cplInterface.readData(meshName3D, dataReadVec, {&wid3, 1}, maxDt, vOut3);

    Vector3d exp0(0.0, 0.0, 10.0);
    Vector3d exp1(0.0, 0.0, 10.0);
    Vector3d exp2(0.0, 0.0, 20.0);
    Vector3d exp3(0.0, 0.0, 20.0);

    BOOST_TEST(vOut0 == exp0);
    BOOST_TEST(vOut1 == exp1);
    BOOST_TEST(vOut2 == exp2);
    BOOST_TEST(vOut3 == exp3);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName3D, dataWriteSca, {vids3D.data(), 4}, pInit);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName3D, dataReadVec, {&wid0, 1}, maxDt, vOut0);
      cplInterface.readData(meshName3D, dataReadVec, {&wid1, 1}, maxDt, vOut1);
      cplInterface.readData(meshName3D, dataReadVec, {&wid2, 1}, maxDt, vOut2);
      cplInterface.readData(meshName3D, dataReadVec, {&wid3, 1}, maxDt, vOut3);

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
