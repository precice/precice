#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "logging/LogMacros.hpp"
#include "mapping/AxialGeoMultiscaleMapping.hpp"
#include "mapping/Mapping.hpp"
#include "math/constants.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(AxialGeoMultiscaleMapping)

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadX)
{
  PRECICE_TEST();
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, coupled along the x-axis.
      Then, the data is mapped from the single vertex of the 1D mesh to defined vertices on the circular inlet of the 3D mesh (hence, "spread").
      The defined vertices are at certain distances from the center, which enables to predict the expected behavior for Hagen-Poiseuille flow.
      Finally, this expected behavior is tested.
  */
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID())); // Point a (1D)
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point A (3D): center, equal to incoming mesh node
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0)); // Point B (3D): distance of 1.0 = r to center
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // Point C (3D): distance of 0.5 = r/2 to center
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X, radius);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(3);
  inValues << 2.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(9);
  outValues = Eigen::VectorXd::Zero(9);

  // Map data
  mapping.computeMapping();
  mapping.computeMapping(); // Check (only in this case) if calling computeMapping() additional times works.
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);
  mapping.computeMapping(); // Check (only in this case) if calling computeMapping() additional times works.

  // Point A (3D): Check if x axis data is doubled at center node (parabolic profile)
  BOOST_TEST(outValues(0) == 2 * inSample.values(0));
  BOOST_TEST(outValues(1) == 0);
  BOOST_TEST(outValues(2) == 0);

  // Point B (3D): Check if x axis data at distance = r is equal to zero
  BOOST_TEST(outValues(3) == 0.0);
  BOOST_TEST(outValues(4) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);

  // Point C (3D): Check if x axis data at distance = r/2 is 3/2 times invalue data
  BOOST_TEST(outValues(6) == 1.5 * inSample.values(0));
  BOOST_TEST(outValues(7) == 0.0);
  BOOST_TEST(outValues(8) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadParabolicZ)
{
  PRECICE_TEST();
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, coupled along the z-axis.
      Then, the data is mapped from the single vertex of the 1D mesh to defined vertices on the circular inlet of the 3D mesh (hence, "spread").
      The defined vertices are at certain distances from the center, which enables to predict the expected behavior for Hagen-Poiseuille flow.
      Finally, this expected behavior is tested.
  */
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point a (1D)
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point A (3D): center, equal to incoming mesh node
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0)); // Point B (3D): distance of 1.0 = r to center
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // Point C (3D): distance of 0.5 = r/2 to center
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(3);
  inValues << 0.0, 0.0, 2.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(9);
  outValues = Eigen::VectorXd::Zero(9);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Point A (3D): Check if x axis data is doubled at center node
  BOOST_TEST(outValues(0) == 0.0);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 2 * inSample.values(2));

  // Point B (3D): Check if x axis data at distance = r is equal to zero
  BOOST_TEST(outValues(3) == 0.0);
  BOOST_TEST(outValues(4) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);

  // Point C (3D): Check if x axis data at distance = r/2 is 3/2 times invalue data
  BOOST_TEST(outValues(6) == 0.0);
  BOOST_TEST(outValues(7) == 0.0);
  BOOST_TEST(outValues(8) == 1.5 * inSample.values(2));
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadUniformZ)
{
  PRECICE_TEST();
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, coupled along the z-axis.
      Then, the data is mapped from the single vertex of the 1D mesh to defined vertices on the circular inlet of the 3D mesh (hence, "spread").
      The values are spread following a uniform spread, it is, all vertices from the 3D mesh receive the same values.
      Finally, this expected behavior is tested.
  */
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point a (1D)
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point A (3D): center, equal to incoming mesh node
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0)); // Point B (3D): distance of 1.0 = r to center
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // Point C (3D): distance of 0.5 = r/2 to center
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(3);
  inValues << 0.0, 0.0, 4.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(9);
  outValues = Eigen::VectorXd::Zero(9);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Point A (3D): Check if z axis data is 4.0
  BOOST_TEST(outValues(0) == 0.0);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 4.0);

  // Point B (3D): Check if z axis data is 4.0
  BOOST_TEST(outValues(3) == 0.0);
  BOOST_TEST(outValues(4) == 0.0);
  BOOST_TEST(outValues(5) == 4.0);

  // Point C (3D): Check if z axis data is 4.0
  BOOST_TEST(outValues(6) == 0.0);
  BOOST_TEST(outValues(7) == 0.0);
  BOOST_TEST(outValues(8) == 4.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectX)
{
  PRECICE_TEST();
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, coupled along the x-axis.
      Then, the data is mapped from multiple defined vertices on the circular inlet of the 3D mesh to the single vertex of the 1D mesh (hence, "collect").
      The defined vertices are at certain distances from the center, which enables to predict the expected behavior for Hagen-Poiseuille flow.
      Finally, this expected behavior is tested.
  */
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point a (3D): center
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0)); // Point b (3D): distance of 1.0 = r to center
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // Point c (3D): distance of 0.5 = r/2 to center
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point A (1D): equal to center of incoming mesh (averaging)
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X, radius);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(9);
  inValues << 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(3);
  outValues = Eigen::VectorXd::Zero(3);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Point A (1D): Check if data is averaged at center node
  BOOST_TEST(outValues(0) == (1 / 3.0) * (inSample.values(0) + inSample.values(3) + inSample.values(6)));
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectZ)
{
  PRECICE_TEST();
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, coupled along the z-axis.
      Then, the data is mapped from multiple defined vertices on the circular inlet of the 3D mesh to the single vertex of the 1D mesh (hence, "collect").
      The defined vertices are at certain distances from the center, which enables to predict the expected behavior for Hagen-Poiseuille flow.
      Finally, this expected behavior is tested.
  */
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point a (3D): center
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0)); // Point b (3D): distance of 1.0 = r to center
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // Point c (3D): distance of 0.5 = r/2 to center
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point A (1D): equal to center of incoming mesh
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z, radius);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(9);
  inValues << 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(3);
  outValues = Eigen::VectorXd::Zero(3);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Point A (1D): Check if data is averaged at center node
  BOOST_TEST(outValues(0) == 0.0);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == (1 / 3.0) * (inSample.values(2) + inSample.values(5) + inSample.values(8)));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
