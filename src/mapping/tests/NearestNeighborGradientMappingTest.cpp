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
  PtrData inDataScalar = inMesh->createData("InDataScalar", 1, 0_dataID);
  PtrData inDataVector = inMesh->createData("InDataVector", 2, 1_dataID);
  inDataScalar->requireDataGradient();
  inDataVector->requireDataGradient();
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
  Eigen::MatrixXd &inGradientsScalar = inDataScalar->gradients();
  Eigen::MatrixXd &inGradientsVector = inDataVector->gradients();
  inGradientsScalar.setOnes();
  inGradientsVector.setOnes();

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
  inGradientsScalar.setZero();
  inGradientsVector.setZero();
  outVertex0.setCoords(inVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  outVertex1.setCoords(inVertex1.getCoords() + Eigen::Vector2d::Constant(0.1));

  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));
  mapping.map(inDataVectorID, outDataVectorID);
  Eigen::Vector4d expected(1.0, 2.0, 3.0, 4.0);
  BOOST_CHECK(equals(expected, outValuesVector));

  // Map data with almost coinciding vertices, should be a little different with the gradient optimization
  inGradientsScalar.setOnes();
  inGradientsVector.setOnes();
  outVertex0.setCoords(inVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  outVertex1.setCoords(inVertex1.getCoords() + Eigen::Vector2d::Constant(0.1));

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
  PtrData inDataScalar = inMesh->createData("InDataScalar", 1, 0_dataID);
  PtrData inDataVector = inMesh->createData("InDataVector", 2, 1_dataID);
  inDataScalar->requireDataGradient();
  inDataVector->requireDataGradient();
  int inDataScalarID = inDataScalar->getID();
  int inDataVectorID = inDataVector->getID();
  inMesh->createVertex(Eigen::Vector2d::Constant(0.0));
  inMesh->createVertex(Eigen::Vector2d::Constant(1.0));

  // Create data
  inMesh->allocateDataValues();
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  Eigen::VectorXd &inValuesVector = inDataVector->values();
  inValuesScalar << 1.0, 2.0;
  inValuesVector << 1.0, 2.0, 3.0, 4.0;

  // Create corresponding gradient data (all gradient values = const = 1)
  Eigen::MatrixXd &inGradientsScalar = inDataScalar->gradients();
  Eigen::MatrixXd &inGradientsVector = inDataVector->gradients();

  inGradientsScalar.col(0) << 2.0, 3.0;
  inGradientsScalar.col(1) << 2.0, 3.0;

  inGradientsVector.col(0) << 2.0, 3.0;
  inGradientsVector.col(1) << 4.0, 5.0;
  inGradientsVector.col(2) << 2.0, 3.0;
  inGradientsVector.col(3) << 4.0, 5.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1, 2_dataID);
  PtrData outDataVector   = outMesh->createData("OutDataVector", 2, 3_dataID);
  int     outDataScalarID = outDataScalar->getID();
  int     outDataVectorID = outDataVector->getID();
  outMesh->createVertex(Eigen::Vector2d::Constant(0.1));
  outMesh->createVertex(Eigen::Vector2d::Constant(1.1));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborGradientMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);

  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0) + 0.5);
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1) + 0.5);

  mapping.map(inDataVectorID, outDataVectorID);
  Eigen::Vector4d        expected(1.5, 2.9, 3.5, 4.9);
  const Eigen::VectorXd &outValuesVector = outDataVector->values();
  BOOST_CHECK(equals(expected, outValuesVector));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
