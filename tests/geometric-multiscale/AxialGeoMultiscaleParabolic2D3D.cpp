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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleParabolic2D3D)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  constexpr int DIM = 3;
  constexpr int N2D = 2;
  constexpr int N3D = 4;

  if (context.isNamed("Fluid2D")) {

    // ---------------- 2D SIDE (WRITES SCALAR, READS VECTOR) ----------------
    auto meshName = "Mesh2D";

    std::array<VertexID, N2D> vids2D;

    // 2D vertices along x at z=0: x = 0, 1
    {
      Vector3d pos0(0.0, 0.0, 0.0);
      Vector3d pos1(1.0, 0.0, 0.0);
      vids2D[0] = cplInterface.setMeshVertex(meshName, pos0);
      vids2D[1] = cplInterface.setMeshVertex(meshName, pos1);
    }

    auto dataAName = "PressureLike"; // scalar (spread to 3D)
    auto dataBName = "VelocityLike"; // vector (collected from 3D)

    // Scalars per 2D vertex
    std::array<double, N2D> valueDataA{};
    // Flattened vectors: [x0,y0,z0, x1,y1,z1]
    std::array<double, N2D * DIM> valueDataB{};

    // Provide initial scalar data if required (for initial SPREAD 2D->3D)
    if (cplInterface.requiresInitialData()) {
      // We choose:
      //   2D vertex x=0: 6
      //   2D vertex x=1: 9
      valueDataA[0] = 6.0;
      valueDataA[1] = 9.0;
      cplInterface.writeData(
          meshName,
          dataAName,
          {vids2D.data(), N2D},
          {valueDataA.data(), N2D});
    }

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // First: collect VelocityLike from 3D -> 2D
    cplInterface.readData(
        meshName,
        dataBName,
        {vids2D.data(), N2D},
        maxDt,
        {valueDataB.data(), static_cast<int>(valueDataB.size())});

    {
      // Geometric situation on 3D side (see below):
      //   3D vertices at x=0: z = 0, 0.5   with Uz = 6, 10
      //   3D vertices at x=1: z = 0.5, 1.0 with Uz = 8, 4
      //
      // COLLECT (3D->2D) bands (nearest in x,z):
      //   2D vertex at x=0 collects from 3D: z=0 and z=0.5 → avg(6,10) = 8
      //   2D vertex at x=1 collects from 3D: z=0.5 and z=1 → avg(8,4) = 6
      //
      // So, on 2D side we expect:
      //   2D(x=0) -> (0,0,8)
      //   2D(x=1) -> (0,0,6)

      int ix0 = 0 * DIM;
      int ix1 = 1 * DIM;

      BOOST_TEST(valueDataB[ix0 + 0] == 0.0);
      BOOST_TEST(valueDataB[ix0 + 1] == 0.0);
      BOOST_TEST(valueDataB[ix0 + 2] == 8.0);

      BOOST_TEST(valueDataB[ix1 + 0] == 0.0);
      BOOST_TEST(valueDataB[ix1 + 1] == 0.0);
      BOOST_TEST(valueDataB[ix1 + 2] == 6.0);
    }

    while (cplInterface.isCouplingOngoing()) {

      // Keep writing the same scalar values on 2D (to be SPREAD to 3D):
      //   2D vertex at x=0: 6.0
      //   2D vertex at x=1: 9.0
      valueDataA[0] = 6.0;
      valueDataA[1] = 9.0;

      cplInterface.writeData(
          meshName,
          dataAName,
          {vids2D.data(), N2D},
          {valueDataA.data(), N2D});

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      // Read collected VelocityLike again: should remain (0,0,8) and (0,0,6)
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
      BOOST_TEST(valueDataB[ix0 + 2] == 8.0);

      BOOST_TEST(valueDataB[ix1 + 0] == 0.0);
      BOOST_TEST(valueDataB[ix1 + 1] == 0.0);
      BOOST_TEST(valueDataB[ix1 + 2] == 6.0);
    }

    cplInterface.finalize();

  } else {

    // ---------------- 3D SIDE (WRITES VECTOR, READS SCALAR) ----------------
    BOOST_TEST(context.isNamed("Fluid3D"));

    auto meshName = "Mesh3D";

    std::array<VertexID, N3D> vids3D;

    {
      // 3D vertices:
      //   x=0: z = 0, 0.5
      //   x=1: z = 0.5, 1.0
      Vector3d pos0(0.0, 0.0, 0.0);
      Vector3d pos1(0.0, 0.0, 0.5);
      Vector3d pos2(1.0, 0.0, 0.5);
      Vector3d pos3(1.0, 0.0, 1.0);
      vids3D[0] = cplInterface.setMeshVertex(meshName, pos0);
      vids3D[1] = cplInterface.setMeshVertex(meshName, pos1);
      vids3D[2] = cplInterface.setMeshVertex(meshName, pos2);
      vids3D[3] = cplInterface.setMeshVertex(meshName, pos3);
    }

    auto dataAName = "PressureLike"; // scalar (spread from 2D)
    auto dataBName = "VelocityLike"; // vector written here

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    // Initial 3D vectors (only z-component non-zero):
    //   x=0: Uz = 6, 10
    //   x=1: Uz = 8, 4
    std::array<double, N3D * DIM> valueDataB3D = {
        0.0, 0.0, 6.0,  // 3D v0
        0.0, 0.0, 10.0, // 3D v1
        0.0, 0.0, 8.0,  // 3D v2
        0.0, 0.0, 4.0,  // 3D v3
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

    // After SPREAD (2D -> 3D) of 2D scalars [6, 9], with a parabolic 2D-3D
    // model and geometry:
    //
    // For x=0 (inputs at z=0, outputs z=0 and 0.5):
    //   distances: r = 0, 0.5 → R_0 = 0.5
    //   r_hat: 0, 1
    //   u(z=0)   = (4/3)*6*(1 - 0)   = 8
    //   u(z=0.5) = (4/3)*6*(1 - 1)   = 0
    //
    // For x=1 (inputs at z=0, outputs z=0.5, 1):
    //   distances: r = 0.5, 1.0 → R_1 = 1
    //   r_hat: 0.5, 1
    //   u(z=0.5) = (4/3)*9*(1 - 0.25) = 9
    //   u(z=1.0) = 0
    //
    // => Expected scalar on 3D:
    //    [8, 0, 9, 0]

    BOOST_TEST(valueDataA3D[0] == 8.0);
    BOOST_TEST(valueDataA3D[1] == 0.0);
    BOOST_TEST(valueDataA3D[2] == 9.0);
    BOOST_TEST(valueDataA3D[3] == 0.0);

    while (cplInterface.isCouplingOngoing()) {

      // Keep writing same 3D vectors so collected 2D result remains:
      //   2D(x=0) -> (0,0,8)
      //   2D(x=1) -> (0,0,6)
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
          0.0,
          0.0,
          4.0,
      };

      cplInterface.writeData(
          meshName,
          dataBName,
          {vids3D.data(), N3D},
          {valueDataB3D.data(), static_cast<int>(valueDataB3D.size())});

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(
          meshName,
          dataAName,
          {vids3D.data(), N3D},
          maxDt,
          {valueDataA3D.data(), N3D});

      BOOST_TEST(valueDataA3D[0] == 8.0);
      BOOST_TEST(valueDataA3D[1] == 0.0);
      BOOST_TEST(valueDataA3D[2] == 9.0);
      BOOST_TEST(valueDataA3D[3] == 0.0);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // GeometricMultiscale
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
