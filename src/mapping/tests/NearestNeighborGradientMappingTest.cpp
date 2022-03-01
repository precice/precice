#include <Eigen/Core>
#include <algorithm>
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborGradientMapping.hpp"
#include "math/constants.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(NearestNeighborGradientMapping)

BOOST_AUTO_TEST_CASE(ConsistentNonIncremental)
{
  PRECICE_TEST(1_rank)
  int dimensions = 2;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createDataWithGradient("InDataScalar", 1, dimensions, 0_dataID);
  PtrData inDataVector   = inMesh->createDataWithGradient("InDataVector", 2, dimensions, 1_dataID);
  int     inDataScalarID = inDataScalar->getID();
  int     inDataVectorID = inDataVector->getID();
  Vertex &inVertex0      = inMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &inVertex1      = inMesh->createVertex(Eigen::Vector2d::Constant(1.0));

  // Create data
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  Eigen::VectorXd &inValuesVector = inDataVector->values();
  inValuesScalar << 1.0, 2.0;
  inValuesVector << 1.0, 2.0, 3.0, 4.0;

  // Create corresponding gradient data (all gradient values = const = 1)
  Eigen::MatrixXd &inGradValuesScalar = inDataScalar->gradientValues();
  Eigen::MatrixXd &inGradValuesVector = inDataVector->gradientValues();
  inGradValuesScalar.setOnes();
  inGradValuesVector.setOnes();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1, 2_dataID);
  PtrData outDataVector   = outMesh->createData("OutDataVector", 2, 3_dataID);
  int     outDataScalarID = outDataScalar->getID();
  int     outDataVectorID = outDataVector->getID();
  Vertex &outVertex0      = outMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &outVertex1      = outMesh->createVertex(Eigen::Vector2d::Constant(1.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborGradientMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  // Distance between in and out vertices is zero
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));

  mapping.map(inDataVectorID, outDataVectorID);
  const Eigen::VectorXd &outValuesVector = outDataVector->values();
  BOOST_CHECK(equals(inValuesVector, outValuesVector));

  // Map data with almost coinciding vertices, with a null gradient, has to result in equal values
  inGradValuesScalar.setZero();
  inGradValuesVector.setZero();
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  inVertex1.setCoords(outVertex1.getCoords() + Eigen::Vector2d::Constant(0.1));

  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));
  mapping.map(inDataVectorID, outDataVectorID);
  Eigen::Vector4d expected(1.0, 2.0, 3.0, 4.0);
  BOOST_CHECK(equals(expected, outValuesVector));

  // Map data with almost coinciding vertices, should be a little different with the gradient optimization
  inGradValuesScalar.setOnes();
  inGradValuesVector.setOnes();
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  inVertex1.setCoords(outVertex1.getCoords() + Eigen::Vector2d::Constant(0.1));

  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0) + 0.2);
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1) + 0.2);
  mapping.map(inDataVectorID, outDataVectorID);
  expected << 1.2, 2.2, 3.2, 4.2;
  BOOST_CHECK(equals(expected, outValuesVector));
}

BOOST_AUTO_TEST_CASE(ConsistentGradientNotConstant)
{

  PRECICE_TEST(1_rank)
  int dimensions = 2;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createDataWithGradient("InDataScalar", 1, dimensions, 0_dataID);
  PtrData inDataVector   = inMesh->createDataWithGradient("InDataVector", 2, dimensions, 1_dataID);
  int     inDataScalarID = inDataScalar->getID();
  int     inDataVectorID = inDataVector->getID();
  Vertex &inVertex0      = inMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &inVertex1      = inMesh->createVertex(Eigen::Vector2d::Constant(1.0));

  // Create data
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  Eigen::VectorXd &inValuesVector = inDataVector->values();
  inValuesScalar << 1.0, 2.0;
  inValuesVector << 1.0, 2.0, 3.0, 4.0;

  // Create corresponding gradient data (all gradient values = const = 1)
  Eigen::MatrixXd &inGradValuesScalar = inDataScalar->gradientValues();
  Eigen::MatrixXd &inGradValuesVector = inDataVector->gradientValues();

  inGradValuesScalar.col(0) << 2.0, 3.0;
  inGradValuesScalar.col(1) << 2.0, 3.0;

  inGradValuesVector.col(0) << 2.0, 3.0;
  inGradValuesVector.col(1) << 5.0, 6.0;
  inGradValuesVector.col(2) << 2.0, 3.0;
  inGradValuesVector.col(3) << 5.0, 6.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1, 2_dataID);
  PtrData outDataVector   = outMesh->createData("OutDataVector", 2, 3_dataID);
  int     outDataScalarID = outDataScalar->getID();
  int     outDataVectorID = outDataVector->getID();
  Vertex &outVertex0      = outMesh->createVertex(Eigen::Vector2d::Constant(0.1));
  Vertex &outVertex1      = outMesh->createVertex(Eigen::Vector2d::Constant(1.1));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborGradientMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0) - 0.5);
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1) - 0.5);

  mapping.map(inDataVectorID, outDataVectorID);
  Eigen::Vector4d        expected(0.5, 0.9, 2.5, 2.9);
  const Eigen::VectorXd &outValuesVector = outDataVector->values();
  BOOST_CHECK(equals(expected, outValuesVector));
}

BOOST_AUTO_TEST_CASE(ConservativeNonIncremental)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inData    = inMesh->createDataWithGradient("InData", 1, dimensions, 0_dataID);
  int     inDataID  = inData->getID();
  Vertex &inVertex0 = inMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &inVertex1 = inMesh->createVertex(Eigen::Vector2d::Constant(1.0));

  // Create data
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValues = inData->values();
  inValues << 1.0, 2.0;

  // Create corresponding gradient
  Eigen::MatrixXd &inGradValues = inData->gradientValues();
  inGradValues.setOnes();

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outData    = outMesh->createData("OutData", 1, 1_dataID);
  int     outDataID  = outData->getID();
  Vertex &outVertex0 = outMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  Vertex &outVertex1 = outMesh->createVertex(Eigen::Vector2d::Constant(1.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborGradientMapping mapping(mapping::Mapping::CONSERVATIVE, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  Eigen::VectorXd &outValues = outData->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == inValues(0));
  BOOST_TEST(outValues(1) == inValues(1));
  outValues = Eigen::VectorXd::Constant(outValues.size(), 0.0);

  // Map data with almost coinciding vertices, has to result in a slight difference that's the value of the gradient optimization
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  inVertex1.setCoords(outVertex1.getCoords() + Eigen::Vector2d::Constant(0.1));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == inValues(0) + 0.2);
  BOOST_TEST(outValues(1) == inValues(1) + 0.2);
  outValues = Eigen::VectorXd::Constant(outValues.size(), 0.0);

  // Map data with coinciding output vertices, has to result in double values with gradient optimization.
  inVertex0.setCoords(Eigen::Vector2d::Constant(0.0));
  inVertex1.setCoords(Eigen::Vector2d::Constant(1.0));

  outVertex1.setCoords(Eigen::Vector2d::Constant(-1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValues(0) == inValues(0) + inValues(1) + 2); //dist*Grad(0,0 - 0,0) = 0 + dist*Grad(0,0 - 0,1) = (1 1) * (1 1)^T
  BOOST_TEST(outValues(1) == 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
