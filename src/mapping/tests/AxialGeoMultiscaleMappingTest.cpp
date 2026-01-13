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

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
