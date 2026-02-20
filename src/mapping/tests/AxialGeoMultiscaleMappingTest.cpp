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
BOOST_AUTO_TEST_CASE(ConsistentSpreadVectorParabolicZ1D3D)
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
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC, mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);
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
BOOST_AUTO_TEST_CASE(ConsistentSpreadVectorUniformX1D3D)
{
  PRECICE_TEST();
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, coupled along the x-axis.
      Then, the data is mapped from the single vertex of the 1D mesh to defined vertices on the circular inlet of the 3D mesh (hence, "spread").
      The values are spread following a uniform spread, meaning that all vertices from the 3D mesh receive the same values.
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
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0)); // Point B (3D): distance of 1.0 = r to center
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // Point C (3D): distance of 0.5 = r/2 to center
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM, mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(3);
  inValues << 4.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(9);
  outValues = Eigen::VectorXd::Zero(9);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Check that the 3D data is initialized to zero and the MultiscaleAxis-component is assigned the same value on all points
  // Point A (3D)
  BOOST_TEST(outValues(0) == 4.0);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);

  // Point B (3D)
  BOOST_TEST(outValues(3) == 4.0);
  BOOST_TEST(outValues(4) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);

  // Point C (3D)
  BOOST_TEST(outValues(6) == 4.0);
  BOOST_TEST(outValues(7) == 0.0);
  BOOST_TEST(outValues(8) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectVectorX1D3D)
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
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X, radius);
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
BOOST_AUTO_TEST_CASE(ConsistentCollectZ1D3D)
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
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z, radius);
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

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarParabolicY1D3D)
{
  PRECICE_TEST();
  // 1D -> 3D, scalar data, parabolic spread along the Y axis
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point a (1D)
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point A (3D): center, equal to incoming mesh node
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0)); // Point B (3D): distance of 1.0 = r to center
  outMesh->createVertex(Eigen::Vector3d(0.5, 0.0, 0.0)); // Point C (3D): distance of 0.5 = r/2 to center
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC, mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(1);
  inValues << 30.0;
  const time::Sample inSample{1, inValues};
  Eigen::VectorXd    outValues(3);
  outValues = Eigen::VectorXd::Zero(3);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // Point A (3D): Check if point data is doubled at center node
  BOOST_TEST(outValues(0) == 2 * inSample.values(0));
  // Point B (3D): Check if point data at distance = r is equal to zero
  BOOST_TEST(outValues(1) == 0);
  // Point C (3D): Check if point data at distance = r/2 is 3/2 times invalue data
  BOOST_TEST(outValues(2) == 1.5 * inSample.values(0));
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarUniformX1D3D)
{
  PRECICE_TEST();
  // 1D -> 3D, scalar data, uniform spread along the X axis
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point a (1D)
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // Point A (3D): center, equal to incoming mesh node
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0)); // Point B (3D): distance of 1.0 = r to center
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // Point C (3D): distance of 0.5 = r/2 to center
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM, mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(1);
  inValues << 42.0;
  const time::Sample inSample{1, inValues};
  Eigen::VectorXd    outValues(3);
  outValues = Eigen::VectorXd::Zero(3);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // All outputs equal to the scalar input
  BOOST_TEST(outValues(0) == 42.0);
  BOOST_TEST(outValues(1) == 42.0);
  BOOST_TEST(outValues(2) == 42.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectScalarZ1D3D)
{
  PRECICE_TEST();
  // 3D -> 1D, scalar data: average over input vertices to the single output vertex, along Z axis
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
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z, radius);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(3);
  inValues << 3.0, 7.0, 5.0;
  const time::Sample inSample{1, inValues};
  Eigen::VectorXd    outValues(1);
  outValues = Eigen::VectorXd::Zero(1);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  // All outputs equal to the scalar input
  BOOST_TEST(outValues(0) == (3.0 + 7.0 + 5.0) / 3);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarUniformZ2D3D)
{
  PRECICE_TEST();
  // 2D -> 3D, scalar data, uniform spread along the Z axis
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, -1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, -1.1, 0.0)); // w0 -> nearest v0
  outMesh->createVertex(Eigen::Vector3d(0.0, -0.2, 0.0)); // w1 -> nearest v1
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.3, 0.0));  // w2 -> nearest v1 (banding)
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.9, 0.0));  // w3 -> nearest v2
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM, mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(3);
  inValues << 10.0, 20.0, 30.0;
  const time::Sample inSample{1, inValues};
  Eigen::VectorXd    outValues(4);
  outValues = Eigen::VectorXd::Zero(4);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outValues(0) == 10.0);
  BOOST_TEST(outValues(1) == 20.0);
  BOOST_TEST(outValues(2) == 20.0);
  BOOST_TEST(outValues(3) == 30.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarParabolicZ2D3D)
{
  PRECICE_TEST();
  // 2D -> 3D, scalar data, uniform spread along the Z axis
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.5, 0.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC, mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(2);
  inValues << 10.0, 20.0;
  const time::Sample inSample{1, inValues};
  Eigen::VectorXd    outValues(4);
  outValues = Eigen::VectorXd::Zero(4);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outValues(0) == 40.0 / 3.0);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 20.0);
  BOOST_TEST(outValues(3) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadUniformVectorY2D3D)
{
  PRECICE_TEST();
  // 2D -> 3D, scalar data, uniform spread along the Y axis
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(-1.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(-1.1, 0.0, 0.0)); // w0 -> nearest v0
  outMesh->createVertex(Eigen::Vector3d(-0.2, 0.0, 0.0)); // w1 -> nearest v1
  outMesh->createVertex(Eigen::Vector3d(0.3, 0.0, 0.0));  // w2 -> nearest v1 (banding)
  outMesh->createVertex(Eigen::Vector3d(0.9, 0.0, 0.0));  // w3 -> nearest v2
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM, mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(9);
  inValues << 0.0, 10.0, 0.0,
      0.0, 20.0, 0.0,
      0.0, 30.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(12);
  outValues = Eigen::VectorXd::Zero(12);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outValues(1) == 10.0);
  BOOST_TEST(outValues(4) == 20.0);
  BOOST_TEST(outValues(7) == 20.0);
  BOOST_TEST(outValues(10) == 30.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadParabolicVectorY2D3D)
{
  PRECICE_TEST();
  // 2D -> 3D, scalar data, uniform spread along the Y axis
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.5));
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.5));
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y, radius, mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC, mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(6);
  inValues << 0.0, 10.0, 0.0,
      0.0, 20.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(12);
  outValues = Eigen::VectorXd::Zero(12);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outValues(1) == 40.0 / 3.0);
  BOOST_TEST(outValues(4) == 0.0);
  BOOST_TEST(outValues(7) == 20.0);
  BOOST_TEST(outValues(10) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectScalarX2D3D)
{
  PRECICE_TEST();
  // 3D -> 2D, scalar data, collect along the X axis
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -1.1));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.9));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.2));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.8));
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -1.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X, radius);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(4);
  inValues << 10.0, 20.0, 30.0, 40.0;
  const time::Sample inSample{1, inValues};
  Eigen::VectorXd    outValues(3);
  outValues = Eigen::VectorXd::Zero(3);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outValues(0) == 15.0);
  BOOST_TEST(outValues(1) == 30.0);
  BOOST_TEST(outValues(2) == 40.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectVectorY2D3D)
{
  PRECICE_TEST();
  // 3D -> 2D, scalar data, collect along the X axis
  constexpr int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -1.1));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.9));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.2));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.8));
  inMesh->allocateDataValues();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -1.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  outMesh->allocateDataValues();

  double radius = 1.0; // radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3, mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT, mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y, radius);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Create data to map
  Eigen::VectorXd inValues(12);
  inValues << 0.0, 10.0, 0.0,
      0.0, 20.0, 0.0,
      0.0, 30.0, 0.0,
      0.0, 40.0, 0.0;
  const time::Sample inSample{3, inValues};
  Eigen::VectorXd    outValues(9);
  outValues = Eigen::VectorXd::Zero(9);

  // Map data
  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outValues(1) == 15.0);
  BOOST_TEST(outValues(4) == 30.0);
  BOOST_TEST(outValues(7) == 40.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarParabolicZ1D2D)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  // 1D mesh (single vertex)
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();

  // 2D mesh (we place points such that distance-from-center is in Y)
  // center: y=0, wall: y=h, mid: y=h/2
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0)); // center
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0)); // wall (h)
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // mid (h/2)
  outMesh->allocateDataValues();

  double radius = 1.0; // interpret as half-gap h for parallel plates

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D2,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(1);
  inValues << 30.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(3);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // 2D Poiseuille between parallel plates:
  // center: (3/2)*avg, wall: 0, mid (h/2): (9/8)*avg
  BOOST_TEST(outValues(0) == 1.5 * inSample.values(0));
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 1.125 * inSample.values(0));
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarUniformZ1D2D)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  outMesh->allocateDataValues();

  double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D2,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(1);
  inValues << 42.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(3);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(0) == 42.0);
  BOOST_TEST(outValues(1) == 42.0);
  BOOST_TEST(outValues(2) == 42.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadVectorParabolicY1D2D)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0)); // center
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0)); // wall
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // mid
  outMesh->allocateDataValues();

  double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D2,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(3);
  inValues << 0.0, 8.0, 0.0; // only Y component
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(9);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Check that only Y-components are populated according to 2D Poiseuille factors
  // Vertex 0 (center)
  BOOST_TEST(outValues(0) == 0.0);
  BOOST_TEST(outValues(1) == 1.5 * inSample.values(1));
  BOOST_TEST(outValues(2) == 0.0);

  // Vertex 1 (wall)
  BOOST_TEST(outValues(3) == 0.0);
  BOOST_TEST(outValues(4) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);

  // Vertex 2 (mid)
  BOOST_TEST(outValues(6) == 0.0);
  BOOST_TEST(outValues(7) == 1.125 * inSample.values(1));
  BOOST_TEST(outValues(8) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadVectorUniformY1D2D)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  outMesh->allocateDataValues();

  double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D2,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::CIRCLE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(3);
  inValues << 0.0, 3.0, 0.0;
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(9);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // All points get the same vector
  BOOST_TEST(outValues(1) == 3.0);
  BOOST_TEST(outValues(4) == 3.0);
  BOOST_TEST(outValues(7) == 3.0);

  // other components remain 0
  BOOST_TEST(outValues(0) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);
  BOOST_TEST(outValues(3) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);
  BOOST_TEST(outValues(6) == 0.0);
  BOOST_TEST(outValues(8) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectScalarX1D2D)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  outMesh->allocateDataValues();

  double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D2,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(3);
  inValues << 3.0, 7.0, 5.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(1);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(0) == (3.0 + 7.0 + 5.0) / 3.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectVectorX1D2D)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  outMesh->allocateDataValues();

  double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D2,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(9);
  inValues << 10.0, 0.0, 0.0,
      20.0, 0.0, 0.0,
      30.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(3);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(0) == (10.0 + 20.0 + 30.0) / 3.0);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarParabolicZ1D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;
  const auto    tol        = boost::test_tools::tolerance(1e-12);

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();

  // Axis=Z => transverse coords are X,Y. Include a corner to test product term.
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0)); // center (s1=0,s2=0)
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0)); // edge   (s1=1,s2=0)
  outMesh->createVertex(Eigen::Vector3d(0.5, 0.0, 0.0)); // mid    (s1=0.5,s2=0)
  outMesh->createVertex(Eigen::Vector3d(0.5, 0.5, 0.0)); // corner (s1=0.5,s2=0.5)
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(1);
  inValues << 30.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(4);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  constexpr double factor = 2.096;
  constexpr double m      = 0.879;

  const double u = inSample.values(0);

  const double b_center = 1.0;             // 1 - 0^2
  const double b_edge   = 0.0;             // 1 - 1^2
  const double b_mid    = 1.0 - 0.5 * 0.5; // 0.75
  const double b_corner = 1.0 - 0.5 * 0.5; // 0.75

  const double expected_center = factor * u * std::pow(b_center, m) * std::pow(b_center, m);
  const double expected_edge   = factor * u * std::pow(b_edge, m) * std::pow(b_center, m); // = 0
  const double expected_mid    = factor * u * std::pow(b_mid, m) * std::pow(b_center, m);
  const double expected_corner = factor * u * std::pow(b_corner, m) * std::pow(b_corner, m);

  BOOST_TEST(outValues(0) == expected_center, tol);
  BOOST_TEST(outValues(1) == 0.0, tol);
  BOOST_TEST(outValues(2) == expected_mid, tol);
  BOOST_TEST(outValues(3) == expected_corner, tol);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarUniformZ1D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.5, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.5, 0.5, 0.0));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(1);
  inValues << 42.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(4);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(0) == 42.0);
  BOOST_TEST(outValues(1) == 42.0);
  BOOST_TEST(outValues(2) == 42.0);
  BOOST_TEST(outValues(3) == 42.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadVectorParabolicX1D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;
  const auto    tol        = boost::test_tools::tolerance(1e-12);

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();

  // Axis=X => transverse coords are Y,Z. Include a corner point (0.5,0.5) in YZ.
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0)); // center
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0)); // edge (s1=1)
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // mid  (s1=0.5)
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.5)); // corner (s1=0.5,s2=0.5)
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(3);
  inValues << 2.0, 0.0, 0.0; // only X component
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(12);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  constexpr double factor = 2.096;
  constexpr double m      = 0.879;

  const double u = inSample.values(0);

  const double b_center = 1.0;
  const double b_edge   = 0.0;
  const double b_mid    = 1.0 - 0.5 * 0.5; // 0.75
  const double b_corner = 1.0 - 0.5 * 0.5; // 0.75

  const double expected_center = factor * u * std::pow(b_center, m) * std::pow(b_center, m);
  const double expected_edge   = 0.0;
  const double expected_mid    = factor * u * std::pow(b_mid, m) * std::pow(b_center, m);
  const double expected_corner = factor * u * std::pow(b_corner, m) * std::pow(b_corner, m);

  // vertex 0
  BOOST_TEST(outValues(0) == expected_center, tol);
  BOOST_TEST(outValues(1) == 0.0, tol);
  BOOST_TEST(outValues(2) == 0.0, tol);
  // vertex 1
  BOOST_TEST(outValues(3) == expected_edge, tol);
  BOOST_TEST(outValues(4) == 0.0, tol);
  BOOST_TEST(outValues(5) == 0.0, tol);
  // vertex 2
  BOOST_TEST(outValues(6) == expected_mid, tol);
  BOOST_TEST(outValues(7) == 0.0, tol);
  BOOST_TEST(outValues(8) == 0.0, tol);
  // vertex 3
  BOOST_TEST(outValues(9) == expected_corner, tol);
  BOOST_TEST(outValues(10) == 0.0, tol);
  BOOST_TEST(outValues(11) == 0.0, tol);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadVectorUniformX1D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.5));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(3);
  inValues << 4.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(12);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  for (int v = 0; v < 4; ++v) {
    BOOST_TEST(outValues(3 * v + 0) == 4.0);
    BOOST_TEST(outValues(3 * v + 1) == 0.0);
    BOOST_TEST(outValues(3 * v + 2) == 0.0);
  }
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectScalarUniformX1D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.5));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(4);
  inValues << 1.0, 2.0, 3.0, 4.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(1);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(0) == (1.0 + 2.0 + 3.0 + 4.0) / 4.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectVectorUniformX1D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.5));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(12);
  inValues << 1.0, 0.0, 0.0,
      2.0, 0.0, 0.0,
      3.0, 0.0, 0.0,
      4.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(3);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(0) == (1.0 + 2.0 + 3.0 + 4.0) / 4.0);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectScalarParabolicX1D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.5));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(4);
  inValues << 10.0, 20.0, 30.0, 40.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(1);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(0) == (10.0 + 20.0 + 30.0 + 40.0) / 4.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectVectorParabolicX1D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.5));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D1D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(12);
  inValues << 10.0, 0.0, 0.0,
      20.0, 0.0, 0.0,
      30.0, 0.0, 0.0,
      40.0, 0.0, 0.0;
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(3);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(0) == (10.0 + 20.0 + 30.0 + 40.0) / 4.0);
  BOOST_TEST(outValues(1) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarUniformZ2D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, -1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, -1.1, 0.0)); // -> v0
  outMesh->createVertex(Eigen::Vector3d(0.0, -0.2, 0.0)); // -> v1
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.3, 0.0));  // -> v1
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.9, 0.0));  // -> v2
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(3);
  inValues << 10.0, 20.0, 30.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(4);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(0) == 10.0);
  BOOST_TEST(outValues(1) == 20.0);
  BOOST_TEST(outValues(2) == 20.0);
  BOOST_TEST(outValues(3) == 30.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadScalarParabolicZ2D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;
  const auto    tol        = boost::test_tools::tolerance(1e-12);

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  // line along Y -> _lineCoord = 1
  inMesh->createVertex(Eigen::Vector3d(0.0, -0.5, 0.0)); // v0
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));  // v1
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  // near v0
  outMesh->createVertex(Eigen::Vector3d(0.0, -0.5, 0.0)); // dist=0
  outMesh->createVertex(Eigen::Vector3d(0.5, -0.5, 0.0)); // dist=0.5
  // near v1
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // dist=0
  outMesh->createVertex(Eigen::Vector3d(0.5, 0.5, 0.0)); // dist=0.5
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Z,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(2);
  inValues << 10.0, 20.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(4);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Expected from your D2D3 square parabolic formula
  constexpr double umax_over_umean = 2.096;
  constexpr double line_factor     = 1.5;
  constexpr double m               = 0.879;
  constexpr double eps             = 1e-12;

  // For v0 and v1: s1 = y_in / radius = +/- 0.5 -> b1raw = 1 - 0.25 = 0.75
  const double s1    = 0.5;
  const double b1raw = 1.0 - s1 * s1; // 0.75
  const double b1    = std::max(eps, b1raw);

  const double pref = (umax_over_umean / line_factor) * std::pow(b1, m - 1.0);

  // dist=0 -> s2=0 -> b2=1
  const double b2_0 = 1.0;
  // dist=0.5 -> s2=0.5 -> b2=0.75
  const double s2_h = 0.5;
  const double b2_h = std::max(0.0, 1.0 - s2_h * s2_h); // 0.75

  const double expected_v0_d0 = inSample.values(0) * pref * std::pow(b2_0, m);
  const double expected_v0_dh = inSample.values(0) * pref * std::pow(b2_h, m);
  const double expected_v1_d0 = inSample.values(1) * pref * std::pow(b2_0, m);
  const double expected_v1_dh = inSample.values(1) * pref * std::pow(b2_h, m);

  BOOST_TEST(outValues(0) == expected_v0_d0, tol);
  BOOST_TEST(outValues(1) == expected_v0_dh, tol);
  BOOST_TEST(outValues(2) == expected_v1_d0, tol);
  BOOST_TEST(outValues(3) == expected_v1_dh, tol);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadVectorUniformY2D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(-1.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(-1.1, 0.0, 0.0)); // -> v0
  outMesh->createVertex(Eigen::Vector3d(-0.2, 0.0, 0.0)); // -> v1
  outMesh->createVertex(Eigen::Vector3d(0.3, 0.0, 0.0));  // -> v1
  outMesh->createVertex(Eigen::Vector3d(0.9, 0.0, 0.0));  // -> v2
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(9);
  inValues << 0.0, 10.0, 0.0,
      0.0, 20.0, 0.0,
      0.0, 30.0, 0.0;
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(12);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  BOOST_TEST(outValues(1) == 10.0);
  BOOST_TEST(outValues(4) == 20.0);
  BOOST_TEST(outValues(7) == 20.0);
  BOOST_TEST(outValues(10) == 30.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentSpreadVectorParabolicY2D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;
  const auto    tol        = boost::test_tools::tolerance(1e-12);

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, -0.5, 0.0)); // v0
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));  // v1
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, -0.5, 0.0)); // dist=0
  outMesh->createVertex(Eigen::Vector3d(0.5, -0.5, 0.0)); // dist=0.5
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));  // dist=0
  outMesh->createVertex(Eigen::Vector3d(0.5, 0.5, 0.0));  // dist=0.5
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::SPREAD,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(6);
  inValues << 0.0, 10.0, 0.0,
      0.0, 20.0, 0.0;
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(12);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  constexpr double umax_over_umean = 2.096;
  constexpr double line_factor     = 1.5;
  constexpr double m               = 0.879;
  constexpr double eps             = 1e-12;

  const double s1    = 0.5;
  const double b1raw = 1.0 - s1 * s1; // 0.75
  const double b1    = std::max(eps, b1raw);
  const double pref  = (umax_over_umean / line_factor) * std::pow(b1, m - 1.0);

  const double b2_0 = 1.0;
  const double s2_h = 0.5;
  const double b2_h = std::max(0.0, 1.0 - s2_h * s2_h); // 0.75

  const double expected_v0_d0 = 10.0 * pref * std::pow(b2_0, m);
  const double expected_v0_dh = 10.0 * pref * std::pow(b2_h, m);
  const double expected_v1_d0 = 20.0 * pref * std::pow(b2_0, m);
  const double expected_v1_dh = 20.0 * pref * std::pow(b2_h, m);

  // Check only Y components: indices 1,4,7,10
  BOOST_TEST(outValues(1) == expected_v0_d0, tol);
  BOOST_TEST(outValues(4) == expected_v0_dh, tol);
  BOOST_TEST(outValues(7) == expected_v1_d0, tol);
  BOOST_TEST(outValues(10) == expected_v1_dh, tol);

  // Other components should remain zero
  BOOST_TEST(outValues(0) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);
  BOOST_TEST(outValues(3) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);
  BOOST_TEST(outValues(6) == 0.0);
  BOOST_TEST(outValues(8) == 0.0);
  BOOST_TEST(outValues(9) == 0.0);
  BOOST_TEST(outValues(11) == 0.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectScalarUniformX2D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.6));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.4));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.2));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.6));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  // span along Z -> _lineCoord = 2
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.5));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.5));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(4);
  inValues << 10.0, 20.0, 30.0, 40.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(3);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Bands: [-0.6,-0.4] -> -0.5 => avg 15; [0.2] -> 0 => 30; [0.6] -> 0.5 => 40
  BOOST_TEST(outValues(0) == 15.0);
  BOOST_TEST(outValues(1) == 30.0);
  BOOST_TEST(outValues(2) == 40.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectVectorUniformY2D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.6));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.4));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.2));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.6));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.5));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.5));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::UNIFORM,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(12);
  inValues << 0.0, 10.0, 0.0,
      0.0, 20.0, 0.0,
      0.0, 30.0, 0.0,
      0.0, 40.0, 0.0;
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(9);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Y-components at outputs: 15, 30, 40
  BOOST_TEST(outValues(1) == 15.0);
  BOOST_TEST(outValues(4) == 30.0);
  BOOST_TEST(outValues(7) == 40.0);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectScalarParabolicX2D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;
  const auto    tol        = boost::test_tools::tolerance(1e-12);

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.6));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.4));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.2));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.6));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  // Keep within [-1,1] so b1 > 0
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.5));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.5));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::X,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(4);
  inValues << 10.0, 20.0, 30.0, 40.0;
  const time::Sample inSample{1, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(3);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Base averages (same as uniform case)
  const double avg0 = 15.0;
  const double avg1 = 30.0;
  const double avg2 = 40.0;

  // Square parabolic collect scaling: 1.03745 * pow(b1, 0.121), b1=max(0,1-s1^2), s1=z/radius
  constexpr double scaleC = 1.03745;
  constexpr double expP   = 0.121;

  const double z0 = -0.5, z1 = 0.0, z2 = 0.5;

  const double b10 = std::max(0.0, 1.0 - z0 * z0); // 0.75
  const double b11 = std::max(0.0, 1.0 - z1 * z1); // 1.0
  const double b12 = std::max(0.0, 1.0 - z2 * z2); // 0.75

  const double expected0 = avg0 * scaleC * std::pow(b10, expP);
  const double expected1 = avg1 * scaleC * std::pow(b11, expP);
  const double expected2 = avg2 * scaleC * std::pow(b12, expP);

  BOOST_TEST(outValues(0) == expected0, tol);
  BOOST_TEST(outValues(1) == expected1, tol);
  BOOST_TEST(outValues(2) == expected2, tol);
}

PRECICE_TEST_SETUP(1_rank);
BOOST_AUTO_TEST_CASE(ConsistentCollectVectorParabolicY2D3D_Square)
{
  PRECICE_TEST();
  constexpr int dimensions = 3;
  const auto    tol        = boost::test_tools::tolerance(1e-12);

  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.6));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.4));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.2));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.6));
  inMesh->allocateDataValues();

  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, -0.5));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.5));
  outMesh->allocateDataValues();

  const double radius = 1.0;

  precice::mapping::AxialGeoMultiscaleMapping mapping(
      mapping::Mapping::CONSISTENT, dimensions,
      mapping::AxialGeoMultiscaleMapping::MultiscaleDimension::D2D3,
      mapping::AxialGeoMultiscaleMapping::MultiscaleType::COLLECT,
      mapping::AxialGeoMultiscaleMapping::MultiscaleAxis::Y,
      radius,
      mapping::AxialGeoMultiscaleMapping::SpreadProfile::PARABOLIC,
      mapping::AxialGeoMultiscaleMapping::MultiscaleCrossSection::SQUARE);

  mapping.setMeshes(inMesh, outMesh);

  Eigen::VectorXd inValues(12);
  inValues << 0.0, 10.0, 0.0,
      0.0, 20.0, 0.0,
      0.0, 30.0, 0.0,
      0.0, 40.0, 0.0;
  const time::Sample inSample{3, inValues};

  Eigen::VectorXd outValues = Eigen::VectorXd::Zero(9);

  mapping.computeMapping();
  mapping.map(inSample, outValues);

  // Base averages for Y component
  const double avg0 = 15.0;
  const double avg1 = 30.0;
  const double avg2 = 40.0;

  constexpr double scaleC = 1.03745;
  constexpr double expP   = 0.121;

  const double z0 = -0.5, z1 = 0.0, z2 = 0.5;

  const double b10 = std::max(0.0, 1.0 - z0 * z0); // 0.75
  const double b11 = std::max(0.0, 1.0 - z1 * z1); // 1.0
  const double b12 = std::max(0.0, 1.0 - z2 * z2); // 0.75

  const double expected0 = avg0 * scaleC * std::pow(b10, expP);
  const double expected1 = avg1 * scaleC * std::pow(b11, expP);
  const double expected2 = avg2 * scaleC * std::pow(b12, expP);

  // Check only Y-components
  BOOST_TEST(outValues(1) == expected0, tol);
  BOOST_TEST(outValues(4) == expected1, tol);
  BOOST_TEST(outValues(7) == expected2, tol);

  // other components remain 0
  BOOST_TEST(outValues(0) == 0.0);
  BOOST_TEST(outValues(2) == 0.0);
  BOOST_TEST(outValues(3) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);
  BOOST_TEST(outValues(6) == 0.0);
  BOOST_TEST(outValues(8) == 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
