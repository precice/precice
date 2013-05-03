// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "RadialBasisFctMappingTest.hpp"
#include "mapping/RadialBasisFctMapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Data.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include "tarch/la/ScalarOperations.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mapping::tests::RadialBasisFctMappingTest)

namespace precice {
namespace mapping {
namespace tests {

tarch::logging::Log RadialBasisFctMappingTest::
  _log ( "precice::mapping::tests::RadialBasisFctMappingTest" );

RadialBasisFctMappingTest:: RadialBasisFctMappingTest()
:
  TestCase ( "precice::mapping::tests::RadialBasisFctMappingTest" )
{}

void RadialBasisFctMappingTest:: run()
{
  testMethod(testThinPlateSplines);
  testMethod(testMultiquadrics);
  testMethod(testInverseMultiquadrics);
  testMethod(testVolumeSplines);
  testMethod(testGaussian);
  testMethod(testCompactThinPlateSplinesC2);
  testMethod(testCompactPolynomialC0);
  testMethod(testCompactPolynomialC6);
}

void RadialBasisFctMappingTest:: testThinPlateSplines()
{
  preciceTrace ( "testThinPlateSplines" );
  ThinPlateSplines fct;
  RadialBasisFctMapping<ThinPlateSplines> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<ThinPlateSplines> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<ThinPlateSplines> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<ThinPlateSplines> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void RadialBasisFctMappingTest:: testMultiquadrics()
{
  preciceTrace ( "testMultiquadrics" );
  Multiquadrics fct(1e-3);
  RadialBasisFctMapping<Multiquadrics> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<Multiquadrics> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<Multiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<Multiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void RadialBasisFctMappingTest:: testInverseMultiquadrics()
{
  preciceTrace ( "testInverseMultiquadrics" );
  InverseMultiquadrics fct(1e-3);
  RadialBasisFctMapping<InverseMultiquadrics> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<InverseMultiquadrics> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<InverseMultiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<InverseMultiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void RadialBasisFctMappingTest:: testVolumeSplines()
{
  preciceTrace ( "testVolumeSplines" );
  VolumeSplines fct;
  RadialBasisFctMapping<VolumeSplines> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<VolumeSplines> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<VolumeSplines> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<VolumeSplines> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void RadialBasisFctMappingTest:: testGaussian()
{
  preciceTrace ( "testGaussian" );
  Gaussian fct(1.0);
  RadialBasisFctMapping<Gaussian> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  RadialBasisFctMapping<Gaussian> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  RadialBasisFctMapping<Gaussian> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  RadialBasisFctMapping<Gaussian> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void RadialBasisFctMappingTest:: testCompactThinPlateSplinesC2()
{
  preciceTrace ( "testCompactThinPlateSplinesC2" );
  double supportRadius = 1.2;
  CompactThinPlateSplinesC2 fct(supportRadius);
  typedef RadialBasisFctMapping<CompactThinPlateSplinesC2> Mapping;
  Mapping consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void RadialBasisFctMappingTest:: testCompactPolynomialC0()
{
  preciceTrace ( "testCompactPolynomialC0" );
  double supportRadius = 1.2;
  CompactPolynomialC0 fct(supportRadius);
  typedef RadialBasisFctMapping<CompactPolynomialC0> Mapping;
  Mapping consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void RadialBasisFctMappingTest:: testCompactPolynomialC6()
{
  preciceTrace ( "testCompactPolynomialC6" );
  double supportRadius = 1.2;
  CompactPolynomialC6 fct(supportRadius);
  typedef RadialBasisFctMapping<CompactPolynomialC6> Mapping;
  Mapping consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void RadialBasisFctMappingTest:: perform2DTestConsistentMapping
(
  Mapping& mapping )
{
  preciceTrace ( "perform2DTestConsistentMapping()" );
  int dimensions = 2;
  using utils::Vector2D;

  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  inMesh->createVertex ( Vector2D(0.0, 0.0) );
  inMesh->createVertex ( Vector2D(1.0, 0.0) );
  inMesh->createVertex ( Vector2D(1.0, 1.0) );
  inMesh->createVertex ( Vector2D(0.0, 1.0) );
  inMesh->allocateDataValues ();
  tarch::la::Vector<4,double> assignValues;
  assignList(assignValues) = 1.0, 2.0, 2.0, 1.0;
  utils::DynVector& values = inData->values();
  values = assignValues;

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();
  mesh::Vertex& vertex = outMesh->createVertex ( Vector2D(0.0) );
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes ( inMesh, outMesh );
  validateEquals ( mapping.hasComputedMapping(), false );

  vertex.setCoords ( Vector2D(0.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  double value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.0 );

  vertex.setCoords ( Vector2D(0.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.0 );

  vertex.setCoords ( Vector2D(0.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.0 );

  vertex.setCoords ( Vector2D(1.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 2.0 );

  vertex.setCoords ( Vector2D(1.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 2.0 );

  vertex.setCoords ( Vector2D(1.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 2.0 );

  vertex.setCoords ( Vector2D(0.5, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.5 );

  vertex.setCoords ( Vector2D(0.5, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.5 );

  vertex.setCoords ( Vector2D(0.5, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.5 );
}

void RadialBasisFctMappingTest:: perform2DTestConservativeMapping
(
  Mapping& mapping )
{
  preciceTrace ( "perform2DTestConservativeMapping()" );
  int dimensions = 2;
  using utils::Vector2D;

  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  mesh::Vertex& vertex0 = inMesh->createVertex ( Vector2D(0.0) );
  mesh::Vertex& vertex1 = inMesh->createVertex ( Vector2D(0.0) );
  inMesh->allocateDataValues ();
  assignList(inData->values()) = 1.0, 2.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID ();
  outMesh->createVertex ( Vector2D(0.0, 0.0) );
  outMesh->createVertex ( Vector2D(1.0, 0.0) );
  outMesh->createVertex ( Vector2D(1.0, 1.0) );
  outMesh->createVertex ( Vector2D(0.0, 1.0) );
  outMesh->allocateDataValues ();
  utils::DynVector& values = outData->values();

  mapping.setMeshes ( inMesh, outMesh );
  validateEquals ( mapping.hasComputedMapping(), false );

  vertex0.setCoords ( Vector2D(0.5, 0.0) );
  vertex1.setCoords ( Vector2D(0.5, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(values, tarch::la::Vector<4,double>(0.5, 0.5, 1.0, 1.0)) );

  vertex0.setCoords ( Vector2D(0.0, 0.5) );
  vertex1.setCoords ( Vector2D(1.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(values, tarch::la::Vector<4,double>(0.5, 1.0, 1.0, 0.5)) );

  vertex0.setCoords ( Vector2D(0.0, 1.0) );
  vertex1.setCoords ( Vector2D(1.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(values, tarch::la::Vector<4,double>(0.0, 2.0, 0.0, 1.0)) );

  vertex0.setCoords ( Vector2D(0.0, 0.0) );
  vertex1.setCoords ( Vector2D(1.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(values, tarch::la::Vector<4,double>(1.0, 0.0, 2.0, 0.0)) );

  vertex0.setCoords ( Vector2D(0.4, 0.5) );
  vertex1.setCoords ( Vector2D(0.6, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( sum(values), 3.0 );
}

void RadialBasisFctMappingTest:: perform3DTestConsistentMapping
(
  Mapping& mapping )
{
  preciceTrace ( "perform3DTestConsistentMapping()" );
  int dimensions = 3;
  using utils::Vector3D;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false));
  mesh::PtrData inData = inMesh->createData("InData", 1);
  int inDataID = inData->getID();
  inMesh->createVertex(Vector3D(0.0, 0.0, 0.0));
  inMesh->createVertex(Vector3D(1.0, 0.0, 0.0));
  inMesh->createVertex(Vector3D(0.0, 1.0, 0.0));
  inMesh->createVertex(Vector3D(1.0, 1.0, 0.0));
  inMesh->createVertex(Vector3D(0.0, 0.0, 1.0));
  inMesh->createVertex(Vector3D(1.0, 0.0, 1.0));
  inMesh->createVertex(Vector3D(0.0, 1.0, 1.0));
  inMesh->createVertex(Vector3D(1.0, 1.0, 1.0));
  inMesh->allocateDataValues();
  utils::DynVector& values = inData->values();
  assignList(values) = 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false));
  mesh::PtrData outData = outMesh->createData("OutData", 1);
  int outDataID = outData->getID();
  mesh::Vertex& vertex = outMesh->createVertex(Vector3D(0.0));
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  validateEquals(mapping.hasComputedMapping(), false);

  vertex.setCoords(Vector3D(0.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  double value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.0);

  vertex.setCoords(Vector3D(0.0, 0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.0);

  vertex.setCoords(Vector3D(0.5, 0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.0);

  vertex.setCoords(Vector3D(1.0, 0.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.0);

  vertex.setCoords(Vector3D(1.0, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.0);

  vertex.setCoords(Vector3D(0.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 2.0);

  vertex.setCoords(Vector3D(1.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 2.0);

  vertex.setCoords(Vector3D(1.0, 1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 2.0);

  vertex.setCoords(Vector3D(0.5, 0.5, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 2.0);

  vertex.setCoords(Vector3D(0.0, 0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.5);

  vertex.setCoords(Vector3D(1.0, 0.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.5);

  vertex.setCoords(Vector3D(0.0, 1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.5);

  vertex.setCoords(Vector3D(1.0, 1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.5);

  vertex.setCoords(Vector3D(0.5, 0.5, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.5);
}

void RadialBasisFctMappingTest:: perform3DTestConservativeMapping
(
  Mapping& mapping )
{
  preciceTrace ( "perform3DTestConservativeMapping()" );
  int dimensions = 3;
  using utils::Vector3D;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false));
  mesh::PtrData inData = inMesh->createData("InData", 1);
  int inDataID = inData->getID();
  mesh::Vertex& vertex0 = inMesh->createVertex(Vector3D(0.0));
  mesh::Vertex& vertex1 = inMesh->createVertex(Vector3D(0.0));
  inMesh->allocateDataValues();
  assignList(inData->values()) = 1.0, 2.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false));
  mesh::PtrData outData = outMesh->createData("OutData", 1);
  int outDataID = outData->getID();
  outMesh->createVertex(Vector3D(0.0, 0.0, 0.0));
  outMesh->createVertex(Vector3D(1.0, 0.0, 0.0));
  outMesh->createVertex(Vector3D(1.0, 1.0, 0.0));
  outMesh->createVertex(Vector3D(0.0, 1.0, 0.0));
  outMesh->createVertex(Vector3D(0.0, 0.0, 1.0));
  outMesh->createVertex(Vector3D(1.0, 0.0, 1.0));
  outMesh->createVertex(Vector3D(1.0, 1.0, 1.0));
  outMesh->createVertex(Vector3D(0.0, 1.0, 1.0));
  outMesh->allocateDataValues();
  utils::DynVector& values = outData->values();
  double expectedSum = sum(inData->values());

  mapping.setMeshes(inMesh, outMesh);
  validateEquals(mapping.hasComputedMapping(), false);

  vertex0.setCoords(Vector3D(0.5, 0.0, 0.0));
  vertex1.setCoords(Vector3D(0.5, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  validateEquals(mapping.hasComputedMapping(), true);
  validateWithParams1(tarch::la::equals(sum(values), expectedSum), values);

//  vertex0.setCoords ( Vector2D(0.0, 0.5) );
//  vertex1.setCoords ( Vector2D(1.0, 0.5) );
//  mapping.computeMapping ();
//  mapping.map ( inDataID, outDataID );
//  validateEquals ( mapping.hasComputedMapping(), true );
//  validate ( equals(values, tarch::la::Vector<4,double>(0.5, 1.0, 1.0, 0.5)) );
//
//  vertex0.setCoords ( Vector2D(0.0, 1.0) );
//  vertex1.setCoords ( Vector2D(1.0, 0.0) );
//  mapping.computeMapping ();
//  mapping.map ( inDataID, outDataID );
//  validateEquals ( mapping.hasComputedMapping(), true );
//  validate ( equals(values, tarch::la::Vector<4,double>(0.0, 2.0, 0.0, 1.0)) );
//
//  vertex0.setCoords ( Vector2D(0.0, 0.0) );
//  vertex1.setCoords ( Vector2D(1.0, 1.0) );
//  mapping.computeMapping ();
//  mapping.map ( inDataID, outDataID );
//  validateEquals ( mapping.hasComputedMapping(), true );
//  validate ( equals(values, tarch::la::Vector<4,double>(1.0, 0.0, 2.0, 0.0)) );
//
//  vertex0.setCoords ( Vector2D(0.4, 0.5) );
//  vertex1.setCoords ( Vector2D(0.6, 0.5) );
//  mapping.computeMapping ();
//  mapping.map ( inDataID, outDataID );
//  validateEquals ( mapping.hasComputedMapping(), true );
//  validateNumericalEquals ( sum(values), 3.0 );
}

}}} // namespace precice, mapping, tests
