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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleParabolic2D3DReverse)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  constexpr int DIM = 3;
  constexpr int N2D = 2;
  constexpr int N3D = 4;

  if (context.isNamed("Fluid2D")) {

    // ---------------- 2D SIDE (WRITES VECTOR, READS SCALAR) ----------------
    auto meshName = "Mesh2D";

    std::array<VertexID, N2D> vids2D;

    // 2D vertices along x at z = 0: x = 0, 1
    {
      Vector3d pos0(0.0, 0.0, 0.0);
      Vector3d pos1(1.0, 0.0, 0.0);
      vids2D[0] = cplInterface.setMeshVertex(meshName, pos0);
      vids2D[1] = cplInterface.setMeshVertex(meshName, pos1);
    }

    auto dataAName = "PressureLike"; // scalar (to be collected from 3D)
    auto dataBName = "VelocityLike"; // vector (to be spread to 3D)

    // Scalars per 2D vertex (read from 3D)
    std::array<double, N2D> valueDataA{};
    // Flattened vectors per 2D vertex (write to 3D)
    std::array<double, N2D * DIM> valueDataB{};

    // Provide initial vector data if required (for initial SPREAD 2D->3D)
    if (cplInterface.requiresInitialData()) {
      // 2D velocity z-components we want to spread:
      //   at x=0: Uz = 10
      //   at x=1: Uz = 20
      valueDataB[0] = 0.0;  // vx0
      valueDataB[1] = 0.0;  // vy0
      valueDataB[2] = 10.0; // vz0

      valueDataB[3] = 0.0;  // vx1
      valueDataB[4] = 0.0;  // vy1
      valueDataB[5] = 20.0; // vz1

      cplInterface.writeData(
          meshName,
          dataBName,
          {vids2D.data(), N2D},
          {valueDataB.data(), static_cast<int>(valueDataB.size())});
    }

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // First: collect scalar PressureLike from 3D -> 2D
    cplInterface.readData(
        meshName,
        dataAName,
        {vids2D.data(), N2D},
        maxDt,
        {valueDataA.data(), N2D});

    {
      // On 3D side we will define pressures:
      //   x=0: P = 2, 4
      //   x=1: P = 6, 8
      //
      // With collect bands:
      //   2D(x=0) collects from (x=0,z=0) and (x=0,z=0.5): avg(2,4) = 3
      //   2D(x=1) collects from (x=1,z=0.5) and (x=1,z=1): avg(6,8) = 7

      BOOST_TEST(valueDataA[0] == 3.0);
      BOOST_TEST(valueDataA[1] == 7.0);
    }

    while (cplInterface.isCouplingOngoing()) {

      // Keep same 2D velocity to be SPREAD:
      //   at x=0: Uz = 10
      //   at x=1: Uz = 20
      valueDataB[0] = 0.0;
      valueDataB[1] = 0.0;
      valueDataB[2] = 10.0;

      valueDataB[3] = 0.0;
      valueDataB[4] = 0.0;
      valueDataB[5] = 20.0;

      cplInterface.writeData(
          meshName,
          dataBName,
          {vids2D.data(), N2D},
          {valueDataB.data(), static_cast<int>(valueDataB.size())});

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      // Read collected PressureLike again: should stay 3 and 7
      cplInterface.readData(
          meshName,
          dataAName,
          {vids2D.data(), N2D},
          maxDt,
          {valueDataA.data(), N2D});

      BOOST_TEST(valueDataA[0] == 3.0);
      BOOST_TEST(valueDataA[1] == 7.0);
    }

    cplInterface.finalize();

  } else {

    // ---------------- 3D SIDE (WRITES SCALAR, READS VECTOR) ----------------
    BOOST_TEST(context.isNamed("Fluid3D"));

    auto meshName = "Mesh3D";

    std::array<VertexID, N3D> vids3D;

    {
      // Geometry as before:
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

    auto dataAName = "PressureLike"; // scalar (written here, collected on 2D)
    auto dataBName = "VelocityLike"; // vector (read here, spread from 2D)

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    // Initial scalar pressures on 3D:
    //   x=0: P = 2, 4
    //   x=1: P = 6, 8
    std::array<double, N3D> valueDataA3D = {2.0, 4.0, 6.0, 8.0};

    cplInterface.writeData(
        meshName,
        dataAName,
        {vids3D.data(), N3D},
        {valueDataA3D.data(), N3D});

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // Read velocity spread from 2D
    std::array<double, N3D * DIM> valueDataB3D{};
    cplInterface.readData(
        meshName,
        dataBName,
        {vids3D.data(), N3D},
        maxDt,
        {valueDataB3D.data(), static_cast<int>(valueDataB3D.size())});

    {
      // Using the same 2D–3D parabolic SPREAD model (factor 4/3) and geometry:
      //
      // For x=0:
      //   2D Uz = 10
      //   distances: r = 0, 0.5 → R_0 = 0.5
      //   r_hat: 0, 1
      //   Uz(0)   = (4/3)*10*(1 - 0)   = 40/3
      //   Uz(0.5) = (4/3)*10*(1 - 1)   = 0
      //
      // For x=1:
      //   2D Uz = 20
      //   distances: r = 0.5, 1 → R_1 = 1
      //   r_hat: 0.5, 1
      //   Uz(0.5) = (4/3)*20*(1 - 0.25) = 20
      //   Uz(1.0) = 0
      //
      // So expected 3D velocities (z-components):
      //   v0.z = 40/3, v1.z = 0, v2.z = 20, v3.z = 0

      // v0
      BOOST_TEST(valueDataB3D[0 * DIM + 0] == 0.0);
      BOOST_TEST(valueDataB3D[0 * DIM + 1] == 0.0);
      BOOST_TEST(valueDataB3D[0 * DIM + 2] == 40.0 / 3.0);

      // v1
      BOOST_TEST(valueDataB3D[1 * DIM + 0] == 0.0);
      BOOST_TEST(valueDataB3D[1 * DIM + 1] == 0.0);
      BOOST_TEST(valueDataB3D[1 * DIM + 2] == 0.0);

      // v2
      BOOST_TEST(valueDataB3D[2 * DIM + 0] == 0.0);
      BOOST_TEST(valueDataB3D[2 * DIM + 1] == 0.0);
      BOOST_TEST(valueDataB3D[2 * DIM + 2] == 20.0);

      // v3
      BOOST_TEST(valueDataB3D[3 * DIM + 0] == 0.0);
      BOOST_TEST(valueDataB3D[3 * DIM + 1] == 0.0);
      BOOST_TEST(valueDataB3D[3 * DIM + 2] == 0.0);
    }

    while (cplInterface.isCouplingOngoing()) {

      // Keep writing the same 3D pressures so the collected 2D scalars remain:
      //   2D x=0 -> 3,  2D x=1 -> 7
      valueDataA3D = {2.0, 4.0, 6.0, 8.0};

      cplInterface.writeData(
          meshName,
          dataAName,
          {vids3D.data(), N3D},
          {valueDataA3D.data(), N3D});

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      cplInterface.readData(
          meshName,
          dataBName,
          {vids3D.data(), N3D},
          maxDt,
          {valueDataB3D.data(), static_cast<int>(valueDataB3D.size())});

      // Velocity field should stay the same pattern:
      BOOST_TEST(valueDataB3D[0 * DIM + 2] == 40.0 / 3.0);
      BOOST_TEST(valueDataB3D[1 * DIM + 2] == 0.0);
      BOOST_TEST(valueDataB3D[2 * DIM + 2] == 20.0);
      BOOST_TEST(valueDataB3D[3 * DIM + 2] == 0.0);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // GeometricMultiscale
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
