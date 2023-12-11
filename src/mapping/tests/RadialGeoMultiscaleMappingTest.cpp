#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/RadialGeoMultiscaleMapping.hpp"
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
BOOST_AUTO_TEST_SUITE(RadialGeoMultiscaleMapping)

BOOST_AUTO_TEST_CASE(testConsistentSpreadX)
{
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, that intersect along the x-axis.
      Then, the data is mapped from the vertices of the 1D mesh to defined vertices of the 3D mesh (hence, SPREAD).
      The defined vertices are close to a 1D vertex - in the sense that their projection onto the x-axis is a nearest neighbor of a 1D vertex -
      and therefore we can predict which value from the 1D vertex should be assigned.
      Finally, this expected behavior is tested.
  */

  PRECICE_TEST(1_rank);
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point a (1D)
  inMesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0)); // Point b (1D)
  inMesh->createVertex(Eigen::Vector3d(6.0, 0.0, 0.0)); // Point c (1D)
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0)); // Point A (3D): closest to Point a
  outMesh->createVertex(Eigen::Vector3d(0.5, 0.5, 0.0)); // Point B (3D): closest to Point a
  outMesh->createVertex(Eigen::Vector3d(4.0, 2.0, 0.0)); // Point C (3D): closest to Point b
  outMesh->createVertex(Eigen::Vector3d(3.0, 1.0, 0.0)); // Point D (3D): closest to Point b
  outMesh->createVertex(Eigen::Vector3d(7.0, 3.0, 0.0)); // Point E (3D): closest to Point c
  outMesh->createVertex(Eigen::Vector3d(6.0, 2.0, 0.0)); // Point F (3D): closest to Point c
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::RadialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::RadialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::RadialGeoMultiscaleMapping::MultiscaleAxis::X);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(9);
  inValues << 1.0, 0.0, 0.0,
      2.0, 0.0, 0.0,
      3.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(18);
  outValues = Eigen::VectorXd::Zero(18);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Check if data is mapped to closest vertex
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Point A (1D) == Point a (3D)
  BOOST_TEST(outValues(0) == inSample.values(0));
  BOOST_TEST(outValues(1) == inSample.values(1));
  BOOST_TEST(outValues(2) == inSample.values(2));

  // Point B (1D) == Point a (3D)
  BOOST_TEST(outValues(3) == inSample.values(0));
  BOOST_TEST(outValues(4) == inSample.values(1));
  BOOST_TEST(outValues(5) == inSample.values(2));

  // Point C (1D) == Point b (3D)
  BOOST_TEST(outValues(6) == inSample.values(3));
  BOOST_TEST(outValues(7) == inSample.values(4));
  BOOST_TEST(outValues(8) == inSample.values(5));

  // Point D (1D) = Point b (3D)
  BOOST_TEST(outValues(9) == inSample.values(3));
  BOOST_TEST(outValues(10) == inSample.values(4));
  BOOST_TEST(outValues(11) == inSample.values(5));

  // Point E (1D) == Point c (3D)
  BOOST_TEST(outValues(12) == inSample.values(6));
  BOOST_TEST(outValues(13) == inSample.values(7));
  BOOST_TEST(outValues(14) == inSample.values(8));

  // Point F (1D) == Point c (3D)
  BOOST_TEST(outValues(15) == inSample.values(6));
  BOOST_TEST(outValues(16) == inSample.values(7));
  BOOST_TEST(outValues(16) == inSample.values(8));
}

BOOST_AUTO_TEST_CASE(testConsistentSpreadZ)
{
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, that intersect along the z-axis.
      Then, the data is mapped from the vertices of the 1D mesh to defined vertices of the 3D mesh (hence, "spread").
      The defined vertices are close to a 1D vertex - in the sense that their projection onto the z-axis is a nearest neighbor of a 1D vertex -
      and therefore we can predict which value from the 1D vertex should be assigned.
      Finally, this expected behavior is tested.
  */

  PRECICE_TEST(1_rank);
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point a (1D)
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 3.0)); // Point b (1D)
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 6.0)); // Point c (1D)
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 1.0)); // Point A (3D): closest to Point a
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.5)); // Point B (3D): closest to Point a
  outMesh->createVertex(Eigen::Vector3d(0.0, 2.0, 4.0)); // Point C (3D): closest to Point b
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 3.0)); // Point D (3D): closest to Point b
  outMesh->createVertex(Eigen::Vector3d(0.0, 3.0, 7.0)); // Point E (3D): closest to Point c
  outMesh->createVertex(Eigen::Vector3d(0.0, 2.0, 6.0)); // Point F (3D): closest to Point c
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::RadialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::RadialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::RadialGeoMultiscaleMapping::MultiscaleAxis::Z);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(9);
  inValues << 1.0, 0.0, 0.0,
      2.0, 0.0, 0.0,
      3.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(18);
  outValues = Eigen::VectorXd::Zero(18);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Check if data is mapped to closest vertex
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Point A (1D) == Point a (3D)
  BOOST_TEST(outValues(0) == inSample.values(0));
  BOOST_TEST(outValues(1) == inSample.values(1));
  BOOST_TEST(outValues(2) == inSample.values(2));

  // Point B (1D) == Point a (3D)
  BOOST_TEST(outValues(3) == inSample.values(0));
  BOOST_TEST(outValues(4) == inSample.values(1));
  BOOST_TEST(outValues(5) == inSample.values(2));

  // Point C (1D) == Point b (3D)
  BOOST_TEST(outValues(6) == inSample.values(3));
  BOOST_TEST(outValues(7) == inSample.values(4));
  BOOST_TEST(outValues(8) == inSample.values(5));

  // Point D (1D) = Point b (3D)
  BOOST_TEST(outValues(9) == inSample.values(3));
  BOOST_TEST(outValues(10) == inSample.values(4));
  BOOST_TEST(outValues(11) == inSample.values(5));

  // Point E (1D) == Point c (3D)
  BOOST_TEST(outValues(12) == inSample.values(6));
  BOOST_TEST(outValues(13) == inSample.values(7));
  BOOST_TEST(outValues(14) == inSample.values(8));

  // Point F (1D) == Point c (3D)
  BOOST_TEST(outValues(15) == inSample.values(6));
  BOOST_TEST(outValues(16) == inSample.values(7));
  BOOST_TEST(outValues(16) == inSample.values(8));
}

BOOST_AUTO_TEST_CASE(testConsistentCollectX)
{
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, that intersect along the x-axis.
      Then, the data is mapped from the vertices of the 3D mesh to defined vertices of the 1D mesh (hence, "collect").
      The defined vertices are close to a 1D vertex - in the sense that their projection onto the x-axis is a nearest neighbor of a 1D vertex -
      and therefore we can predict that the 1D vertex value should be the mean of all close 3D vertices.
      Finally, this expected behavior is tested.
  */

  PRECICE_TEST(1_rank);
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0)); // Point a (3D)
  inMesh->createVertex(Eigen::Vector3d(1.0, 2.0, 0.0)); // Point b (3D)
  inMesh->createVertex(Eigen::Vector3d(2.0, 1.0, 0.0)); // Point c (3D)
  inMesh->createVertex(Eigen::Vector3d(4.0, 1.0, 0.0)); // Point d (3D)
  inMesh->createVertex(Eigen::Vector3d(5.0, 2.0, 0.0)); // Point e (3D)
  inMesh->createVertex(Eigen::Vector3d(7.0, 1.0, 0.0)); // Point f (3D)
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point A (1D): Averaging from points a,b
  outMesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0)); // Point B (1D): Averaging from points c,d
  outMesh->createVertex(Eigen::Vector3d(6.0, 0.0, 0.0)); // Point C (1D): Averaging from points e,f
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::RadialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::RadialGeoMultiscaleMapping::MultiscaleType::COLLECT, mapping::RadialGeoMultiscaleMapping::MultiscaleAxis::X);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(18);
  inValues << 1.0, 0.0, 0.0,
      3.0, 0.0, 0.0,
      5.0, 0.0, 0.0,
      7.0, 0.0, 0.0,
      2.0, 0.0, 0.0,
      4.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(9);
  outValues = Eigen::VectorXd::Zero(9);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Check if data is mapped to closest vertex
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Point A (1D) == average of Points a,b (3D)
  BOOST_TEST(outValues(0) == (inSample.values(0) + inSample.values(3)) / 2);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);

  // Point B (1D) == average of Points c,d (3D)
  BOOST_TEST(outValues(3) == (inSample.values(6) + inSample.values(9)) / 2);
  BOOST_TEST(outValues(4) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);

  // Point C (1D) == average of Points e,f (3D)
  BOOST_TEST(outValues(6) == (inSample.values(12) + inSample.values(15)) / 2);
  BOOST_TEST(outValues(7) == 0.0);
  BOOST_TEST(outValues(8) == 0.0);
}

BOOST_AUTO_TEST_CASE(testConsistentCollectZ)
{
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, that intersect along the x-axis.
      Then, the data is mapped from the vertices of the 3D mesh to defined vertices of the 1D mesh (hence, COLLECT).
      The defined vertices are close to a 1D vertex and therefore we can predict that the 1D vertex value should be the mean of all close 3D vertices.
      Finally, this expected behavior is tested.
  */

  PRECICE_TEST(1_rank);
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 1.0)); // Point a (3D)
  inMesh->createVertex(Eigen::Vector3d(0.0, 2.0, 1.0)); // Point b (3D)
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 2.0)); // Point c (3D)
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 4.0)); // Point d (3D)
  inMesh->createVertex(Eigen::Vector3d(0.0, 2.0, 5.0)); // Point e (3D)
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 7.0)); // Point f (3D)
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point A (1D): Averaging from points a,b
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 3.0)); // Point B (1D): Averaging from points c,d
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 6.0)); // Point C (1D): Averaging from points e,f
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::RadialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::RadialGeoMultiscaleMapping::MultiscaleType::COLLECT, mapping::RadialGeoMultiscaleMapping::MultiscaleAxis::Z);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(18);
  inValues << 1.0, 0.0, 0.0,
      3.0, 0.0, 0.0,
      5.0, 0.0, 0.0,
      7.0, 0.0, 0.0,
      2.0, 0.0, 0.0,
      4.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(9);
  outValues = Eigen::VectorXd::Zero(9);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Check if data is mapped to closest vertex
  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Point A (1D) == average of Points a,b (3D)
  BOOST_TEST(outValues(0) == (inSample.values(0) + inSample.values(3)) / 2);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);

  // Point B (1D) == average of Points c,d (3D)
  BOOST_TEST(outValues(3) == (inSample.values(6) + inSample.values(9)) / 2);
  BOOST_TEST(outValues(4) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);

  // Point C (1D) == average of Points e,f (3D)
  BOOST_TEST(outValues(6) == (inSample.values(12) + inSample.values(15)) / 2);
  BOOST_TEST(outValues(7) == 0.0);
  BOOST_TEST(outValues(8) == 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
