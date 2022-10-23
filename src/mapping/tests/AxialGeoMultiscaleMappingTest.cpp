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

BOOST_AUTO_TEST_CASE(testConsistentSpread)
{
  PRECICE_TEST(1_rank);
  int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inData    = inMesh->createData("InData", 3, 0_dataID);
  int     inDataID  = inData->getID();
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector3d::Constant(0.0));
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValues = inData->values();
  inValues << 0.0, 0.0, 2.0;

  double radius = 1.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outData    = outMesh->createData("OutData", 3, 2_dataID);
  int     outDataID  = outData->getID();
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // center, equal to incoming mesh node
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0)); // distance of 1.0 = r to center
  Vertex &outVertex2 = outMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // distance of 0.5 = r/2 to center
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::SPREAD, radius);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  const Eigen::VectorXd &outValues = outData->values();

  // Check if data is doubled at center node
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == 2 * inValues(0));
  BOOST_TEST(outValues(1) == 2 * inValues(1));
  BOOST_TEST(outValues(2) == 2 * inValues(2));
  // Check if data at distance = r is equal to zero
  BOOST_TEST(outValues(3) == 0.0);
  BOOST_TEST(outValues(4) == 0.0);
  BOOST_TEST(outValues(5) == 0.0);
  // Check if data at distance = r/2 is 3/2 times invalue data
  BOOST_TEST(outValues(6) == 1.5 * inValues(0));
  BOOST_TEST(outValues(7) == 1.5 * inValues(1));
  BOOST_TEST(outValues(8) == 1.5 * inValues(2));
}

BOOST_AUTO_TEST_CASE(testConsistentCollect)
{
  PRECICE_TEST(1_rank);
  int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createData("InDataScalar", 1, 0_dataID);
  PtrData inDataVector   = inMesh->createData("InDataVector", 3, 1_dataID);
  int     inDataScalarID = inDataScalar->getID();
  int     inDataVectorID = inDataVector->getID();
  Vertex &inVertex0      = inMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // center
  Vertex &inVertex1      = inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0)); // distance of 1.0 = r to center (if circle in x-y plane?)
  Vertex &inVertex2      = inMesh->createVertex(Eigen::Vector3d(0.0, 0.5, 0.0)); // distance of 0.5 = r/2 to center
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  Eigen::VectorXd &inValuesVector = inDataVector->values();
  inValuesScalar << 1.0, 2.0, 3.0;
  inValuesVector << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
  double radius = 1.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1, 2_dataID);
  PtrData outDataVector   = outMesh->createData("OutDataVector", 3, 3_dataID);
  int     outDataScalarID = outDataScalar->getID();
  int     outDataVectorID = outDataVector->getID();
  Vertex &outVertex0      = outMesh->createVertex(Eigen::Vector3d::Constant(0.0)); // equal to center of incoming mesh
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions, mapping::AxialGeoMultiscaleMapping::SPREAD, radius);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  mapping.map(inDataVectorID, outDataVectorID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  const Eigen::VectorXd &outValuesVector = outDataVector->values();

  // Check if scalar data is equal and vector data is halved at center node
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesVector(0) == 0.5 * inValuesVector(0));
  BOOST_TEST(outValuesVector(1) == 0.5 * inValuesVector(1));
  BOOST_TEST(outValuesVector(2) == 0.5 * inValuesVector(2));

  // Check if invalue vector data at distance = r/2 is 2/3 times outvalue vector data
  BOOST_TEST(outValuesVector(6) == (2 / 3) * inValuesVector(0));
  BOOST_TEST(outValuesVector(7) == (2 / 3) * inValuesVector(1));
  BOOST_TEST(outValuesVector(8) == (2 / 3) * inValuesVector(2));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()