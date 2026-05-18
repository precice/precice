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
BOOST_AUTO_TEST_CASE(AxialGeoMultiscaleUniform2D3DReverse)
{
  PRECICE_TEST();

  using Eigen::Vector3d;

  Participant cplInterface(context.name, context.config(), context.rank, context.size);

  constexpr int DIM = 3;
  constexpr int N2D = 2;
  constexpr int N3D = 3;

  if (context.isNamed("Fluid2D")) {

    // ---------------- 2D SIDE (READS SCALAR, WRITES VECTOR) ----------------
    auto meshName = "Mesh2D";

    std::array<VertexID, N2D> vids2D;

    {
      Vector3d pos0(0.0, 0.0, 0.0);
      Vector3d pos1(0.5, 0.0, 0.0);
      vids2D[0] = cplInterface.setMeshVertex(meshName, pos0);
      vids2D[1] = cplInterface.setMeshVertex(meshName, pos1);
    }

    auto dataAName = "PressureLike"; // scalar (read here)
    auto dataBName = "VelocityLike"; // vector (written here)

    std::array<double, N2D>       valueDataA{}; // scalar per 2D vertex
    std::array<double, N2D * DIM> valueDataB{}; // flattened vectors

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    // First read: COLLECT (3D -> 2D) of PressureLike
    cplInterface.readData(
        meshName,
        dataAName,
        vids2D,
        maxDt,
        valueDataA);

    // After COLLECT (3D -> 2D) of [6, 8, 8]:
    // - 2D vertex x=0 gets 3D(x=0)         → 6
    // - 2D vertex x=0.5 gets avg(3D x=0.5,0.5)   → avg(8,8) = 8
    BOOST_TEST(valueDataA[0] == 6.0);
    BOOST_TEST(valueDataA[1] == 8.0);

    while (cplInterface.isCouplingOngoing()) {

      // Velocity values on 2D (to be SPREAD to 3D):
      //   2D vertex x=0: (0,0,6)
      //   2D vertex x=0.5: (0,0,9)
      valueDataB = {
          0.0, 0.0, 6.0, // 2D vertex 0
          0.0, 0.0, 9.0, // 2D vertex 1
      };

      cplInterface.writeData(
          meshName,
          dataBName,
          vids2D,
          valueDataB);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      // Read collected PressureLike again (3D -> 2D):
      // expected to remain [6, 8]
      cplInterface.readData(
          meshName,
          dataAName,
          vids2D,
          maxDt,
          valueDataA);

      BOOST_TEST(valueDataA[0] == 6.0);
      BOOST_TEST(valueDataA[1] == 8.0);
    }

    cplInterface.finalize();

  } else {

    // ---------------- 3D SIDE (WRITES SCALAR, READS VECTOR) ----------------
    BOOST_TEST(context.isNamed("Fluid3D"));

    auto meshName = "Mesh3D";

    std::array<VertexID, N3D> vids3D;

    {
      Vector3d pos0(0.0, 0.0, 0.0);
      Vector3d pos1(0.5, 0.2, 0.0);
      Vector3d pos2(0.5, 0.5, 0.0);
      vids3D[0] = cplInterface.setMeshVertex(meshName, pos0);
      vids3D[1] = cplInterface.setMeshVertex(meshName, pos1);
      vids3D[2] = cplInterface.setMeshVertex(meshName, pos2);
    }

    auto dataAName = "PressureLike"; // scalar written here
    auto dataBName = "VelocityLike"; // vector read here

    BOOST_REQUIRE(cplInterface.requiresInitialData());

    std::array<double, N3D>       valueDataA3D = {6.0, 8.0, 8.0};
    std::array<double, N3D * DIM> valueDataB3D{};

    // Initial write of PressureLike (3D -> 2D)
    cplInterface.writeData(
        meshName,
        dataAName,
        vids3D,
        valueDataA3D);

    cplInterface.initialize();
    double maxDt = cplInterface.getMaxTimeStepSize();

    while (cplInterface.isCouplingOngoing()) {

      // Keep writing the same scalar pattern [6, 8, 8]
      cplInterface.writeData(
          meshName,
          dataAName,
          vids3D,
          valueDataA3D);

      cplInterface.advance(maxDt);
      maxDt = cplInterface.getMaxTimeStepSize();

      // Read SPREAD (2D -> 3D) VelocityLike:
      cplInterface.readData(
          meshName,
          dataBName,
          vids3D,
          maxDt,
          valueDataB3D);

      // After SPREAD (2D -> 3D) of [(0,0,6), (0,0,9)]:
      // - 3D vertex x=0 gets 2D(x=0)       → (0,0,6)
      // - 3D vertices x=0.5,0.5 get 2D(x=0.5)    → (0,0,9)
      //   => [ (0,0,6), (0,0,9), (0,0,9) ]
      int ix0 = 0 * DIM;
      int ix1 = 1 * DIM;
      int ix2 = 2 * DIM;

      BOOST_TEST(valueDataB3D[ix0 + 0] == 0.0);
      BOOST_TEST(valueDataB3D[ix0 + 1] == 0.0);
      BOOST_TEST(valueDataB3D[ix0 + 2] == 6.0);

      BOOST_TEST(valueDataB3D[ix1 + 0] == 0.0);
      BOOST_TEST(valueDataB3D[ix1 + 1] == 0.0);
      BOOST_TEST(valueDataB3D[ix1 + 2] == 9.0);

      BOOST_TEST(valueDataB3D[ix2 + 0] == 0.0);
      BOOST_TEST(valueDataB3D[ix2 + 1] == 0.0);
      BOOST_TEST(valueDataB3D[ix2 + 2] == 9.0);
    }

    cplInterface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // GeometricMultiscale
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
