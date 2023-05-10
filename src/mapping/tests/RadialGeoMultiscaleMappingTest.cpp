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

BOOST_AUTO_TEST_CASE(testConsistentSpread)
{
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, that intersect along the x-axis.
      Then, the data is mapped from the vertices of the 1D mesh to defined vertices of the 3D mesh (hence, SPREAD).
      The defined vertices are close to a 1D vertex and therefore we can predict which value from the 1D vertex should be assigned.
      Finally, this expected behavior is tested.  
  */

  PRECICE_TEST(1_rank);
  int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh                  inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData                  inData    = inMesh->createData("InData", 1, 0_dataID);
  int                      inDataID  = inData->getID();
  [[maybe_unused]] Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  [[maybe_unused]] Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  [[maybe_unused]] Vertex &inVertex2 = inMesh->createVertex(Eigen::Vector3d(6.0, 0.0, 0.0));
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValues = inData->values();
  inValues << 1.0, 2.0, 3.0;

  // Create mesh to map to
  PtrMesh                  outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData                  outData    = outMesh->createData("OutData", 1, 2_dataID);
  int                      outDataID  = outData->getID();
  [[maybe_unused]] Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0)); // closest to first 1D vertex
  [[maybe_unused]] Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector3d(4.0, 2.0, 0.0)); // closest to second 1D vertex
  [[maybe_unused]] Vertex &outVertex2 = outMesh->createVertex(Eigen::Vector3d(7.0, 3.0, 0.0)); // closest to third 1D vertex
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::RadialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::RadialGeoMultiscaleMapping::SPREAD, mapping::RadialGeoMultiscaleMapping::X);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  const Eigen::VectorXd &outValues = outData->values();

  // Check if data is mapped to closest vertex
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == inValues(0));
  BOOST_TEST(outValues(1) == inValues(1));
  BOOST_TEST(outValues(2) == inValues(2));
}

BOOST_AUTO_TEST_CASE(testConsistentCollect)
{
  /*  The following test works by creating two dimensionally heterogeneous meshes, namely 1D and 3D, that intersect along the x-axis.
      Then, the data is mapped from the vertices of the 3D mesh to defined vertices of the 1D mesh (hence, COLLECT).
      The defined vertices are close to a 1D vertex and therefore we can predict that the 1D vertex value should be the mean of all close 3D vertices.
      Finally, this expected behavior is tested.  
  */

  PRECICE_TEST(1_rank);
  int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh                  inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData                  inData    = inMesh->createData("InData", 1, 0_dataID);
  int                      inDataID  = inData->getID();
  [[maybe_unused]] Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  [[maybe_unused]] Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector3d(1.0, 2.0, 0.0));
  [[maybe_unused]] Vertex &inVertex2 = inMesh->createVertex(Eigen::Vector3d(2.0, 1.0, 0.0));
  [[maybe_unused]] Vertex &inVertex3 = inMesh->createVertex(Eigen::Vector3d(4.0, 1.0, 0.0));
  [[maybe_unused]] Vertex &inVertex4 = inMesh->createVertex(Eigen::Vector3d(5.0, 2.0, 0.0));
  [[maybe_unused]] Vertex &inVertex5 = inMesh->createVertex(Eigen::Vector3d(7.0, 1.0, 0.0));
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValues = inData->values();
  inValues << 1.0, 3.0, 5.0, 7.0, 2.0, 4.0;

  // Create mesh to map to
  PtrMesh                  outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData                  outData    = outMesh->createData("OutData", 1, 2_dataID);
  int                      outDataID  = outData->getID();
  [[maybe_unused]] Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  [[maybe_unused]] Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  [[maybe_unused]] Vertex &outVertex2 = outMesh->createVertex(Eigen::Vector3d(6.0, 0.0, 0.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::RadialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::RadialGeoMultiscaleMapping::COLLECT, mapping::RadialGeoMultiscaleMapping::X);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  const Eigen::VectorXd &outValues = outData->values();

  // Check if data is mapped to closest vertex
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == (inValues(0) + inValues(1)) / 2);
  BOOST_TEST(outValues(1) == (inValues(2) + inValues(3)) / 2);
  BOOST_TEST(outValues(2) == (inValues(4) + inValues(5)) / 2);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()