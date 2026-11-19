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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleScalarParabolic2D3D_Square)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("Fluid2D")) {
    const auto meshName2D   = "Mesh2D";
    const auto dataWriteSca = "PressureLike";
    const auto dataReadVec  = "VelocityLike";

    // 2D line along Y (embedded in 3D)
    Vector3d v0(0.0, -0.5, 0.0);
    Vector3d v1(0.0, 0.5, 0.0);

    auto vid0 = cplInterface.setMeshVertex(meshName2D, v0);
    auto vid1 = cplInterface.setMeshVertex(meshName2D, v1);

    std::array<int, 2> vids2D = {vid0, vid1};

    Vector3d valueReadVec = Vector3d::Zero();

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // read collected VelocityLike from 3D (we'll initialize it from 3D side so average is known)
    cplInterface.readData(meshName2D, dataReadVec, {&vid0, 1}, maxDt, valueReadVec);
    Vector3d expectedVec(0.0, 0.0, 6.0);
    BOOST_TEST(valueReadVec == expectedVec);

    while (cplInterface.isCouplingOngoing()) {
      // write scalar values on the 2D mesh (these are the "mean" values which are spread to 3D)
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

    // 3D vertices: two bands (nearest to each 2D vertex), each with one point at dist=0 and dist=0.5
    Vector3d w0(0.0, -0.5, 0.0); // nearest v0, dist=0
    Vector3d w1(0.5, -0.5, 0.0); // nearest v0, dist=0.5
    Vector3d w2(0.0, 0.5, 0.0);  // nearest v1, dist=0
    Vector3d w3(0.5, 0.5, 0.0);  // nearest v1, dist=0.5

    auto wid0 = cplInterface.setMeshVertex(meshName3D, w0);
    auto wid1 = cplInterface.setMeshVertex(meshName3D, w1);
    auto wid2 = cplInterface.setMeshVertex(meshName3D, w2);
    auto wid3 = cplInterface.setMeshVertex(meshName3D, w3);

    std::array<int, 4> vids3D = {wid0, wid1, wid2, wid3};

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    // initialize VelocityLike on 3D so that collect to 2D is predictable: avg scalar component = 6 in z
    Eigen::VectorXd vInit(12);
    vInit << 0.0, 0.0, 8.0,
        0.0, 0.0, 4.0,
        0.0, 0.0, 6.0,
        0.0, 0.0, 6.0;

    cplInterface.writeData(meshName3D, dataWriteVec, vids3D, vInit);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // Read parabolically spread scalar on 3D from 2D, SQUARE cross section
    double p0 = 0.0, p1 = 0.0, p2 = 0.0, p3 = 0.0;

    cplInterface.readData(meshName3D, dataReadSca, {&wid0, 1}, maxDt, {&p0, 1});
    cplInterface.readData(meshName3D, dataReadSca, {&wid1, 1}, maxDt, {&p1, 1});
    cplInterface.readData(meshName3D, dataReadSca, {&wid2, 1}, maxDt, {&p2, 1});
    cplInterface.readData(meshName3D, dataReadSca, {&wid3, 1}, maxDt, {&p3, 1});

    constexpr double umax_over_umean = 2.096;
    constexpr double line_factor     = 1.5;
    constexpr double m               = 0.879;
    constexpr double eps             = 1e-12;

    // In spread D2D3 square formula, s1 = input()->vertex(inIdx).getCoords()[_lineCoord] / radius.
    // Here input line is along Y, so _lineCoord = Y, s1 = +/-0.5.
    const double s1    = 0.5;
    const double b1raw = 1.0 - s1 * s1; // 0.75
    const double b1    = std::max(eps, b1raw);

    const double pref = (umax_over_umean / line_factor) * std::pow(b1, m - 1.0);

    // dist=0 -> s2=0 -> b2=1
    const double s2_0 = 0.0;
    const double b2_0 = std::max(0.0, 1.0 - s2_0 * s2_0); // 1

    // dist=0.5 -> s2=0.5 -> b2=0.75
    const double s2_h = 0.5;
    const double b2_h = std::max(0.0, 1.0 - s2_h * s2_h); // 0.75

    // Input scalar values from Fluid2D are p(v0)=10, p(v1)=20
    const double expected_p0 = 10.0 * pref * std::pow(b2_0, m);
    const double expected_p1 = 10.0 * pref * std::pow(b2_h, m);
    const double expected_p2 = 20.0 * pref * std::pow(b2_0, m);
    const double expected_p3 = 20.0 * pref * std::pow(b2_h, m);

    BOOST_TEST(p0 == expected_p0);
    BOOST_TEST(p1 == expected_p1);
    BOOST_TEST(p2 == expected_p2);
    BOOST_TEST(p3 == expected_p3);

    while (cplInterface.isCouplingOngoing()) {
      cplInterface.writeData(meshName3D, dataWriteVec, vids3D, vInit);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName3D, dataReadSca, {&wid0, 1}, maxDt, {&p0, 1});
      cplInterface.readData(meshName3D, dataReadSca, {&wid1, 1}, maxDt, {&p1, 1});
      cplInterface.readData(meshName3D, dataReadSca, {&wid2, 1}, maxDt, {&p2, 1});
      cplInterface.readData(meshName3D, dataReadSca, {&wid3, 1}, maxDt, {&p3, 1});

      BOOST_TEST(p0 == expected_p0);
      BOOST_TEST(p1 == expected_p1);
      BOOST_TEST(p2 == expected_p2);
      BOOST_TEST(p3 == expected_p3);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
