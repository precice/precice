#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <array>
#include <precice/Participant.hpp>
#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(GeometricMultiscale)

PRECICE_TEST_SETUP("Fluid2D"_on(1_rank), "Fluid3D"_on(1_rank))
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleUniform2D3D)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  constexpr int DIM = 3;
  constexpr int N2D = 2;
  constexpr int N3D = 3;

  if (context.isNamed("Fluid2D")) {

    // ---------------- 2D SIDE (WRITES SCALAR, READS VECTOR) ----------------
    auto meshName = "Mesh2D";

    std::array<VertexID, N2D> vids2D;

    // Explicit 2D vertices (spread/collect along x-direction)
    {
      Vector3d pos0(0.0, 0.0, 0.0);
      Vector3d pos1(1.0, 0.0, 0.0);
      vids2D[0] = cplInterface.setMeshVertex(meshName, pos0);
      vids2D[1] = cplInterface.setMeshVertex(meshName, pos1);
    }

    auto dataAName = "PressureLike"; // scalar
    auto dataBName = "VelocityLike"; // vector (3 components)

    // Scalars per 2D vertex
    std::array<double, N2D> valueDataA{};
    // Flattened vectors: [x0,y0,z0, x1,y1,z1]
    std::array<double, N2D * DIM> valueDataB{};

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    cplInterface.readData(
        meshName,
        dataBName,
        {vids2D.data(), N2D},
        maxDt,
        {valueDataB.data(), static_cast<int>(valueDataB.size())});

    {
      // After COLLECT (3D -> 2D):
      // - 2D vertex at x=0 collects from 3D vertex at x=0  → (0,0,6)
      // - 2D vertex at x=1 collects from 3D vertices x=1,2 → avg(10,8) = (0,0,9)

      // 2D vertex 0
      int ix0 = 0 * DIM;
      BOOST_TEST(valueDataB[ix0 + 0] == 0.0);
      BOOST_TEST(valueDataB[ix0 + 1] == 0.0);
      BOOST_TEST(valueDataB[ix0 + 2] == 6.0);
      // 2D vertex 1
      int ix1 = 1 * DIM;
      BOOST_TEST(valueDataB[ix1 + 0] == 0.0);
      BOOST_TEST(valueDataB[ix1 + 1] == 0.0);
      BOOST_TEST(valueDataB[ix1 + 2] == 9.0);
    }

    while (cplInterface.isCouplingOngoing()) {

      // Pressure values on 2D (to be SPREAD to 3D):
      //   2D vertex at x=0: 6.0
      //   2D vertex at x=1: 8.0
      valueDataA[0] = 6.0;
      valueDataA[1] = 8.0;

      cplInterface.writeData(
          meshName,
          dataAName,
          {vids2D.data(), N2D},
          {valueDataA.data(), N2D});

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      // Read collected VelocityLike again:
      // expected to remain:
      //   2D(x=0) -> (0,0,6)
      //   2D(x=1) -> (0,0,9)
      cplInterface.readData(
          meshName,
          dataBName,
          {vids2D.data(), N2D},
          maxDt,
          {valueDataB.data(), static_cast<int>(valueDataB.size())});

      int ix0 = 0 * DIM;
      int ix1 = 1 * DIM;

      BOOST_TEST(valueDataB[ix0 + 0] == 0.0);
      BOOST_TEST(valueDataB[ix0 + 1] == 0.0);
      BOOST_TEST(valueDataB[ix0 + 2] == 6.0);

      BOOST_TEST(valueDataB[ix1 + 0] == 0.0);
      BOOST_TEST(valueDataB[ix1 + 1] == 0.0);
      BOOST_TEST(valueDataB[ix1 + 2] == 9.0);
    }

    cplInterface.finalize();

  } else {

    // ---------------- 3D SIDE (WRITES VECTOR, READS SCALAR) ----------------
    BOOST_TEST(context.isNamed("Fluid3D"));

    auto meshName = "Mesh3D";

    std::array<VertexID, N3D> vids3D;

    {
      Vector3d pos0(0.0, 0.0, 0.0);
      Vector3d pos1(1.0, 0.0, 1.0);
      Vector3d pos2(2.0, 0.0, 2.0);
      vids3D[0] = cplInterface.setMeshVertex(meshName, pos0);
      vids3D[1] = cplInterface.setMeshVertex(meshName, pos1);
      vids3D[2] = cplInterface.setMeshVertex(meshName, pos2);
    }

    auto dataAName = "PressureLike"; // scalar (spread from 2D)
    auto dataBName = "VelocityLike"; // vector written here

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    std::array<double, N3D * DIM> valueDataB3D = {
        0.0, 0.0, 6.0,  // 3D vertex x=0
        0.0, 0.0, 10.0, // 3D vertex x=1
        0.0, 0.0, 8.0,  // 3D vertex x=2
    };

    cplInterface.writeData(
        meshName,
        dataBName,
        {vids3D.data(), N3D},
        {valueDataB3D.data(), static_cast<int>(valueDataB3D.size())});

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    std::array<double, N3D> valueDataA3D{};
    cplInterface.readData(
        meshName,
        dataAName,
        {vids3D.data(), N3D},
        maxDt,
        {valueDataA3D.data(), N3D});

    // After SPREAD (2D -> 3D) of [6, 8]:
    // - 2D vertex x=0 (6) → 3D vertex x=0  → 6
    // - 2D vertex x=1 (8) → 3D vertices x=1 and x=2 → 8, 8
    BOOST_TEST(valueDataA3D[0] == 6.0);
    BOOST_TEST(valueDataA3D[1] == 8.0);
    BOOST_TEST(valueDataA3D[2] == 8.0);

    while (cplInterface.isCouplingOngoing()) {

      // Keep writing the same 3D values (6, 10, 8) so the collected result stays:
      //   2D(x=0) = 6, 2D(x=1) = 9
      valueDataB3D = {
          0.0,
          0.0,
          6.0,
          0.0,
          0.0,
          10.0,
          0.0,
          0.0,
          8.0,
      };

      cplInterface.writeData(
          meshName,
          dataBName,
          {vids3D.data(), N3D},
          {valueDataB3D.data(), static_cast<int>(valueDataB3D.size())});

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      // After each SPREAD step, we still expect:
      //   3D: [6, 8, 8]
      cplInterface.readData(
          meshName,
          dataAName,
          {vids3D.data(), N3D},
          maxDt,
          {valueDataA3D.data(), N3D});

      BOOST_TEST(valueDataA3D[0] == 6.0);
      BOOST_TEST(valueDataA3D[1] == 8.0);
      BOOST_TEST(valueDataA3D[2] == 8.0);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // GeometricMultiscale
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
