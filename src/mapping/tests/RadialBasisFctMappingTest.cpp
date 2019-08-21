#include "testing/Testing.hpp"

#include "mapping/RadialBasisFctMapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Data.hpp"
#include "mesh/Vertex.hpp"
#include "math/math.hpp"

using namespace precice;
using namespace precice::mapping;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(RadialBasisFunctionMapping)

// Forward declarations, see end of file for definitions
void perform2DTestConsistentMapping(Mapping& mapping);
void perform2DTestConservativeMapping(Mapping& mapping);
void perform3DTestConsistentMapping(Mapping& mapping);
void perform3DTestConservativeMapping(Mapping& mapping);

BOOST_AUTO_TEST_CASE(MapThinPlateSplines)
{
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  ThinPlateSplines fct;
  RadialBasisFctMapping<ThinPlateSplines> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<ThinPlateSplines> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<ThinPlateSplines> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<ThinPlateSplines> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapMultiquadrics)
{
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  Multiquadrics fct(1e-3);
  RadialBasisFctMapping<Multiquadrics> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<Multiquadrics> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<Multiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<Multiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapInverseMultiquadrics)
{
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  InverseMultiquadrics fct(1e-3);
  RadialBasisFctMapping<InverseMultiquadrics> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<InverseMultiquadrics> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<InverseMultiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<InverseMultiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapVolumeSplines)
{
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  VolumeSplines fct;
  RadialBasisFctMapping<VolumeSplines> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<VolumeSplines> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<VolumeSplines> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<VolumeSplines> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapGaussian)
{
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  Gaussian fct(1.0);
  RadialBasisFctMapping<Gaussian> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<Gaussian> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<Gaussian> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<Gaussian> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}


BOOST_AUTO_TEST_CASE(MapThinPlaceSplinesC2)
{
  double supportRadius = 1.2;
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  CompactThinPlateSplinesC2 fct(supportRadius);
  using Mapping = RadialBasisFctMapping<CompactThinPlateSplinesC2>;
  Mapping consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapCompactPolynomialC0)
{
  double supportRadius = 1.2;
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  CompactPolynomialC0 fct(supportRadius);
  using Mapping = RadialBasisFctMapping<CompactPolynomialC0>;
  Mapping consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(MapCompactPolynomialC6)
{
  double supportRadius = 1.2;
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  CompactPolynomialC6 fct(supportRadius);
  using Mapping = RadialBasisFctMapping<CompactPolynomialC6>;
  Mapping consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

BOOST_AUTO_TEST_CASE(DeadAxis2D)
{
  int dimensions = 2;
  
  bool xDead = false;
  bool yDead = true;
  bool zDead = false;

  ThinPlateSplines fct;
  RadialBasisFctMapping<ThinPlateSplines> mapping(Mapping::CONSISTENT, dimensions, fct,
                                                  xDead, yDead, zDead);

  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  inMesh->createVertex ( Eigen::Vector2d(0.0, 1.0) );
  inMesh->createVertex ( Eigen::Vector2d(1.0, 1.0) );
  inMesh->createVertex ( Eigen::Vector2d(2.0, 1.0) );
  inMesh->createVertex ( Eigen::Vector2d(3.0, 1.0) );
  inMesh->allocateDataValues ();
  inData->values() << 1.0, 2.0, 2.0, 1.0;
  
  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();
  mesh::Vertex& vertex = outMesh->createVertex ( Eigen::Vector2d::Zero() );
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes ( inMesh, outMesh );
  BOOST_TEST ( mapping.hasComputedMapping() == false );

  vertex.setCoords ( Eigen::Vector2d(0.0, 3.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  double value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 1.0 );
}

BOOST_AUTO_TEST_CASE(DeadAxis3D)
{
  int dimensions = 3;
  using Eigen::Vector3d;

  double supportRadius = 1.2;
  CompactPolynomialC6 fct(supportRadius);
  bool xDead = false;
  bool yDead = true;
  bool zDead = false;
  using Mapping = RadialBasisFctMapping<CompactPolynomialC6>;
  Mapping mapping(Mapping::CONSISTENT, dimensions, fct, xDead, yDead, zDead);

  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  inMesh->createVertex ( Vector3d(0.0, 3.0, 0.0) );
  inMesh->createVertex ( Vector3d(1.0, 3.0, 0.0) );
  inMesh->createVertex ( Vector3d(0.0, 3.0, 1.0) );
  inMesh->createVertex ( Vector3d(1.0, 3.0, 1.0) );
  inMesh->allocateDataValues ();
  inData->values() << 1.0, 2.0, 3.0, 4.0;
  
  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();
  outMesh->createVertex ( Vector3d(0.0, 2.9, 0.0) );
  outMesh->createVertex ( Vector3d(0.8, 2.9, 0.1) );
  outMesh->createVertex ( Vector3d(0.1, 2.9, 0.9) );
  outMesh->createVertex ( Vector3d(1.1, 2.9, 1.1) );
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes ( inMesh, outMesh );
  BOOST_TEST ( mapping.hasComputedMapping() == false );

  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  BOOST_CHECK ( mapping.hasComputedMapping() );

  BOOST_TEST ( outData->values()[0] = 1.0 );
  BOOST_TEST ( outData->values()[1] = 2.0 );
  BOOST_TEST ( outData->values()[2] = 2.9 );
  BOOST_TEST ( outData->values()[3] = 4.3 );
}

void perform2DTestConsistentMapping(Mapping& mapping )
{
  int dimensions = 2;
  
  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  inMesh->createVertex ( Eigen::Vector2d(0.0, 0.0) );
  inMesh->createVertex ( Eigen::Vector2d(1.0, 0.0) );
  inMesh->createVertex ( Eigen::Vector2d(1.0, 1.0) );
  inMesh->createVertex ( Eigen::Vector2d(0.0, 1.0) );
  inMesh->allocateDataValues ();
  inData->values() << 1.0, 2.0, 2.0, 1.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();
  mesh::Vertex& vertex = outMesh->createVertex ( Eigen::Vector2d::Constant(0.0) );
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes ( inMesh, outMesh );
  BOOST_TEST ( mapping.hasComputedMapping() == false );

  vertex.setCoords ( Eigen::Vector2d(0.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  double value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 1.0 );

  vertex.setCoords ( Eigen::Vector2d(0.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 1.0 );

  vertex.setCoords ( Eigen::Vector2d(0.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 1.0 );

  vertex.setCoords ( Eigen::Vector2d(1.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 2.0 );

  vertex.setCoords ( Eigen::Vector2d(1.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 2.0 );

  vertex.setCoords ( Eigen::Vector2d(1.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 2.0 );

  vertex.setCoords ( Eigen::Vector2d(0.5, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 1.5 );

  vertex.setCoords ( Eigen::Vector2d(0.5, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 1.5 );

  vertex.setCoords ( Eigen::Vector2d(0.5, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( value == 1.5 );
}

void perform2DTestConservativeMapping(Mapping& mapping)
{
  using testing::equals;
  int dimensions = 2;
  
  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  mesh::Vertex& vertex0 = inMesh->createVertex ( Eigen::Vector2d::Zero() );
  mesh::Vertex& vertex1 = inMesh->createVertex ( Eigen::Vector2d::Zero() );
  inMesh->allocateDataValues ();
  inData->values() << 1.0, 2.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID ();
  outMesh->createVertex ( Eigen::Vector2d(0.0, 0.0) );
  outMesh->createVertex ( Eigen::Vector2d(1.0, 0.0) );
  outMesh->createVertex ( Eigen::Vector2d(1.0, 1.0) );
  outMesh->createVertex ( Eigen::Vector2d(0.0, 1.0) );
  outMesh->allocateDataValues ();
  Eigen::VectorXd& values = outData->values();

  mapping.setMeshes ( inMesh, outMesh );
  BOOST_TEST ( mapping.hasComputedMapping() == false );

  vertex0.setCoords ( Eigen::Vector2d(0.5, 0.0) );
  vertex1.setCoords ( Eigen::Vector2d(0.5, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_CHECK ( equals(values, Eigen::Vector4d(0.5, 0.5, 1.0, 1.0)) );

  vertex0.setCoords ( Eigen::Vector2d(0.0, 0.5) );
  vertex1.setCoords ( Eigen::Vector2d(1.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_CHECK ( equals( values, Eigen::Vector4d(0.5, 1.0, 1.0, 0.5)) );

  vertex0.setCoords ( Eigen::Vector2d(0.0, 1.0) );
  vertex1.setCoords ( Eigen::Vector2d(1.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_CHECK ( equals(values, Eigen::Vector4d(0.0, 2.0, 0.0, 1.0)) );

  vertex0.setCoords ( Eigen::Vector2d(0.0, 0.0) );
  vertex1.setCoords ( Eigen::Vector2d(1.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_CHECK ( equals(values, Eigen::Vector4d(1.0, 0.0, 2.0, 0.0)) );

  vertex0.setCoords ( Eigen::Vector2d(0.4, 0.5) );
  vertex1.setCoords ( Eigen::Vector2d(0.6, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  BOOST_CHECK ( mapping.hasComputedMapping() );
  BOOST_TEST ( values.sum() == 3.0 );
}

void perform3DTestConsistentMapping(Mapping& mapping )
{
  int dimensions = 3;
  
  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false));
  mesh::PtrData inData = inMesh->createData("InData", 1);
  int inDataID = inData->getID();
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 1.0));
  inMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  inMesh->allocateDataValues();
  inData->values() << 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false));
  mesh::PtrData outData = outMesh->createData("OutData", 1);
  int outDataID = outData->getID();
  mesh::Vertex& vertex = outMesh->createVertex(Eigen::Vector3d::Zero());
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false );

  vertex.setCoords(Eigen::Vector3d(0.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(0.5, 0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.0);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(1.0, 1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(0.5, 0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 2.0);

  vertex.setCoords(Eigen::Vector3d(0.0, 0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(1.0, 0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(0.0, 1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(1.0, 1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.5);

  vertex.setCoords(Eigen::Vector3d(0.5, 0.5, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(value == 1.5);
}

void perform3DTestConservativeMapping(Mapping& mapping)
{
  int dimensions = 3;
  
  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false));
  mesh::PtrData inData = inMesh->createData("InData", 1);
  int inDataID = inData->getID();
  mesh::Vertex& vertex0 = inMesh->createVertex(Eigen::Vector3d::Zero());
  mesh::Vertex& vertex1 = inMesh->createVertex(Eigen::Vector3d::Zero());
  inMesh->allocateDataValues();
  inData->values() << 1.0, 2.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false));
  mesh::PtrData outData = outMesh->createData("OutData", 1);
  int outDataID = outData->getID();
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  outMesh->createVertex(Eigen::Vector3d(0.0, 1.0, 1.0));
  outMesh->allocateDataValues();
  Eigen::VectorXd& values = outData->values();
  double expectedSum = inData->values().sum();

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  vertex0.setCoords(Eigen::Vector3d(0.5, 0.0, 0.0));
  vertex1.setCoords(Eigen::Vector3d(0.5, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_CHECK(mapping.hasComputedMapping() );
  BOOST_TEST(values.sum() == expectedSum);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
