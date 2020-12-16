#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborGradientMapping.hpp"
#include "math/constants.hpp"
#include "mesh/Data.hpp"
#include "mesh/Gradient.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(NearestNeighborGradientMapping)

BOOST_AUTO_TEST_CASE(Scalar2DConsistent)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createData("InDataScalar", 1);
  PtrGradient inGradientScalar = inMesh->createGradient("InDataScalar", 1);
  
  int     inDataScalarID = inDataScalar->getID();
  int     inGradScalarID = inGradientScalar->getID();
  
  Vertex &inVertex0      = inMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &inVertex1      = inMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  inMesh->allocateDataValues();
  inMesh->allocateGradientValues();
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  inValuesScalar << 1.0, 2.0;
  
  Eigen::MatrixXd &inGradientValuesScalar = inGradientScalar->values();
  inGradientValuesScalar << 1.0, 0.0,
                            1.0, 0.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1);
  int     outDataScalarID = outDataScalar->getID();
  Vertex &outVertex0      = outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &outVertex1      = outMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborGradientMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));

  // Map data with almost coinciding vertices, has to result in values + gradient effect
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d::Constant(0.1));
  inVertex1.setCoords(outVertex1.getCoords() - Eigen::Vector2d::Constant(0.1));
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0) + 0.1);
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1) - 0.1);
}

BOOST_AUTO_TEST_CASE(Vector2DConsistent)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  PtrData         inDataVector = inMesh->createData("InDataVector", 2);
  PtrGradient inGradientVector = inMesh->createGradient("InDataVector", 2);
  
  int     inDataVectorID =     inDataVector->getID();
  int     inGradVectorID = inGradientVector->getID();
  
  Vertex &inVertex0      = inMesh->createVertex(Eigen::Vector2d(0.0, 0.0)); // p_0
  Vertex &inVertex1      = inMesh->createVertex(Eigen::Vector2d(1.0, 1.0)); // p_1
  inMesh->allocateDataValues();
  inMesh->allocateGradientValues();
  Eigen::VectorXd &inValuesVector = inDataVector->values();
  inValuesVector << 1.0, 2.0, // (u, v)_0
                    2.0, 4.0; // (u, v)_1
  
  // u = 1+(du/dx)x+(du/dy)y => du/dx = du/dy = 0.5
  // v = 2+(dv/dx)x+(dv/dy)y => dv/dx = dv/dy = 0.5
  Eigen::MatrixXd &inGradientValuesVector = inGradientVector->values();
  inGradientValuesVector.block(0,0,2,2) << 0.5, 0.5, // [ (du/dx)_0 (du/dy)_0 ]
                                           1.0, 1.0; // [ (dv/dx)_0 (dv/dy)_0 ]

  inGradientValuesVector.block(0,2,2,2) << 0.5, 0.5, // [ (du/dx)_1 (du/dy)_1 ]
                                           1.0, 1.0; // [ (dv/dx)_1 (dv/dy)_1 ]
  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  PtrData outDataVector   = outMesh->createData("OutDataVector", 2);
  int     outDataVectorID = outDataVector->getID();
  Vertex &outVertex0      = outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &outVertex1      = outMesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborGradientMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  mapping.computeMapping();
  mapping.map(inDataVectorID, outDataVectorID);
  const Eigen::VectorXd &outValuesVector = outDataVector->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesVector(0) == inValuesVector(0));
  BOOST_TEST(outValuesVector(1) == inValuesVector(1));
  BOOST_TEST(outValuesVector(2) == inValuesVector(2));
  BOOST_TEST(outValuesVector(3) == inValuesVector(3));

  // Map data with almost coinciding vertices, has to result in values + gradient effect.
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d(0.1, 0.1));
  inVertex1.setCoords(outVertex1.getCoords() - Eigen::Vector2d(0.1, 0.1));
  mapping.computeMapping();
  mapping.map(inDataVectorID, outDataVectorID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesVector(0) == inValuesVector(0) + 0.1);
  BOOST_TEST(outValuesVector(1) == inValuesVector(1) + 0.2);
  BOOST_TEST(outValuesVector(2) == inValuesVector(2) - 0.1);
  BOOST_TEST(outValuesVector(3) == inValuesVector(3) - 0.2);
}

BOOST_AUTO_TEST_CASE(Vector2DData3DConsistent)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  PtrData         inDataVector = inMesh->createData("InDataVector", 3);
  PtrGradient inGradientVector = inMesh->createGradient("InDataVector", 3);
  
  int     inDataVectorID =     inDataVector->getID();
  int     inGradVectorID = inGradientVector->getID();
  
  Vertex &inVertex0      = inMesh->createVertex(Eigen::Vector2d(0.0, 0.0)); // p_0
  Vertex &inVertex1      = inMesh->createVertex(Eigen::Vector2d(1.0, 1.0)); // p_1
  inMesh->allocateDataValues();
  inMesh->allocateGradientValues();
  Eigen::VectorXd &inValuesVector = inDataVector->values();
  inValuesVector << 1.0, 2.0, 0.0, // (u, v, w)_0
                    2.0, 4.0, 6.0; // (u, v, w)_1
  
  // u = 1+(du/dx)x+(du/dy)y => du/dx = du/dy = 0.5
  // v = 2+(dv/dx)x+(dv/dy)y => dv/dx = dv/dy = 0.5
  // w = 0+(dw/dx)x+(dw/dy)y => dw/dx = dw/dy = 3.0
  Eigen::MatrixXd &inGradientValuesVector = inGradientVector->values();
  inGradientValuesVector.block(0,0,3,2) << 0.5, 0.5, // [ (du/dx)_0 (du/dy)_0 ]
                                           1.0, 1.0, // [ (dv/dx)_0 (dv/dy)_0 ]
                                           3.0, 3.0; // [ (dw/dx)_0 (dw/dy)_0 ]

  inGradientValuesVector.block(0,2,3,2) << 0.5, 0.5, // [ (du/dx)_1 (du/dy)_1 ]
                                           1.0, 1.0, // [ (dv/dx)_1 (dv/dy)_1 ]
                                           3.0, 3.0; // [ (dw/dx)_0 (dw/dy)_0 ]
  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  PtrData outDataVector   = outMesh->createData("OutDataVector", 3);
  int     outDataVectorID = outDataVector->getID();
  Vertex &outVertex0      = outMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex &outVertex1      = outMesh->createVertex(Eigen::Vector2d(1.0, 1.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborGradientMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  mapping.computeMapping();
  mapping.map(inDataVectorID, outDataVectorID);
  const Eigen::VectorXd &outValuesVector = outDataVector->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesVector(0) == inValuesVector(0));
  BOOST_TEST(outValuesVector(1) == inValuesVector(1));
  BOOST_TEST(outValuesVector(2) == inValuesVector(2));
  BOOST_TEST(outValuesVector(3) == inValuesVector(3));
  BOOST_TEST(outValuesVector(4) == inValuesVector(4));
  BOOST_TEST(outValuesVector(5) == inValuesVector(5));
  

  // Map data with almost coinciding vertices, has to result in values + gradient effect.
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector2d(0.1, 0.1));
  inVertex1.setCoords(outVertex1.getCoords() - Eigen::Vector2d(0.1, 0.1));
  mapping.computeMapping();
  mapping.map(inDataVectorID, outDataVectorID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesVector(0) == inValuesVector(0) + 0.1);
  BOOST_TEST(outValuesVector(1) == inValuesVector(1) + 0.2);
  BOOST_TEST(outValuesVector(2) == inValuesVector(2) + 0.6);
  BOOST_TEST(outValuesVector(3) == inValuesVector(3) - 0.1);
  BOOST_TEST(outValuesVector(4) == inValuesVector(4) - 0.2);
  BOOST_TEST(outValuesVector(5) == inValuesVector(5) - 0.6);
}

BOOST_AUTO_TEST_CASE(Scalar3DConsistent)
{
  PRECICE_TEST(1_rank);
  int dimensions = 3;
  using testing::equals;

  // Create mesh to map from
  PtrMesh inMesh(new Mesh("InMesh", dimensions, false, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createData("InDataScalar", 1);
  PtrGradient inGradientScalar = inMesh->createGradient("InDataScalar", 1);
  
  int     inDataScalarID = inDataScalar->getID();
  int     inGradScalarID = inGradientScalar->getID();
  
  Vertex &inVertex0      = inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &inVertex1      = inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  inMesh->allocateDataValues();
  inMesh->allocateGradientValues();
  Eigen::VectorXd &inValuesScalar = inDataScalar->values();
  inValuesScalar << 1.0, 2.0;
  
  Eigen::MatrixXd &inGradientValuesScalar = inGradientScalar->values();
  inGradientValuesScalar << 1.0, 0.0, 0.0,
                            1.0, 0.0, 0.0;

  // Create mesh to map to
  PtrMesh outMesh(new Mesh("OutMesh", dimensions, false, testing::nextMeshID()));
  PtrData outDataScalar   = outMesh->createData("OutDataScalar", 1);
  int     outDataScalarID = outDataScalar->getID();
  Vertex &outVertex0      = outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &outVertex1      = outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  precice::mapping::NearestNeighborGradientMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  // Map data with coinciding vertices, has to result in equal values.
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  const Eigen::VectorXd &outValuesScalar = outDataScalar->values();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0));
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1));

  // Map data with almost coinciding vertices, has to result in values + gradient effect
  inVertex0.setCoords(outVertex0.getCoords() + Eigen::Vector3d::Constant(0.1));
  inVertex1.setCoords(outVertex1.getCoords() - Eigen::Vector3d::Constant(0.1));
  mapping.computeMapping();
  mapping.map(inDataScalarID, outDataScalarID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  BOOST_TEST(outValuesScalar(0) == inValuesScalar(0) + 0.1);
  BOOST_TEST(outValuesScalar(1) == inValuesScalar(1) - 0.1);
}


BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
