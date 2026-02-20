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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleVectorParabolic1D2D)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("Fluid1D")) {
    const auto meshName = "Mesh1D";
    Vector3d   posOne   = Vector3d::Constant(0.0);
    auto       vid      = cplInterface.setMeshVertex(meshName, posOne);

    const auto dataAName = "VelocityLike";
    const auto dataBName = "PressureLike";

    double valueDataB = 0.0;

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    cplInterface.readData(meshName, dataBName, {&vid, 1}, maxDt, {&valueDataB, 1});
    const double expectedDataB = 6.0;
    BOOST_TEST(valueDataB == expectedDataB);

    while (cplInterface.isCouplingOngoing()) {
      Vector3d valueDataA(0.0, 0.0, 8.0);
      cplInterface.writeData(meshName, dataAName, {&vid, 1}, valueDataA);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(meshName, dataBName, {&vid, 1}, maxDt, {&valueDataB, 1});
      BOOST_TEST(valueDataB == expectedDataB);
    }

    cplInterface.finalize();

  } else {
    BOOST_TEST(context.isNamed("Fluid2D"));

    const auto meshName  = "Mesh2D";
    const auto dataAName = "VelocityLike";
    const auto dataBName = "PressureLike";

    Vector3d p0(0.0, 0.0, 0.0);
    Vector3d p1(0.0, 0.5, 0.0);

    auto vid0 = cplInterface.setMeshVertex(meshName, p0);
    auto vid1 = cplInterface.setMeshVertex(meshName, p1);

    std::array<int, 2> vids = {vid0, vid1};

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    Eigen::VectorXd initValues(2);
    initValues << 8.0, 4.0;
    cplInterface.writeData(meshName, dataBName, {vids.data(), 2}, initValues);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    while (cplInterface.isCouplingOngoing()) {

      cplInterface.writeData(meshName, dataBName, {vids.data(), 2}, initValues);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      Vector3d vOut0 = Vector3d::Zero();
      Vector3d vOut1 = Vector3d::Zero();

      cplInterface.readData(meshName, dataAName, {&vid0, 1}, maxDt, vOut0);
      cplInterface.readData(meshName, dataAName, {&vid1, 1}, maxDt, vOut1);

      constexpr double factor = 1.5;
      Vector3d         expected0(0.0, 0.0, factor * 8.0 * (1.0 - 0.0 * 0.0));
      Vector3d         expected1(0.0, 0.0, factor * 8.0 * (1.0 - 0.5 * 0.5));

      BOOST_TEST(vOut0 == expected0);
      BOOST_TEST(vOut1 == expected1);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
