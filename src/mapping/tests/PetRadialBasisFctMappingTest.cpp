// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "PetRadialBasisFctMappingTest.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Data.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include "tarch/la/ScalarOperations.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mapping::tests::PetRadialBasisFctMappingTest)

namespace precice {
namespace mapping {
namespace tests {

tarch::logging::Log PetRadialBasisFctMappingTest::
  _log ( "precice::mapping::tests::PetRadialBasisFctMappingTest" );

PetRadialBasisFctMappingTest:: PetRadialBasisFctMappingTest()
:
  TestCase ( "precice::mapping::tests::PetRadialBasisFctMappingTest" ),
  tolerance(1e-5)
{}

void PetRadialBasisFctMappingTest:: run()
{
  testMethod(testPetThinPlateSplines);
  testMethod(testPetMultiquadrics);
  testMethod(testPetInverseMultiquadrics);
  testMethod(testPetVolumeSplines);
  testMethod(testPetGaussian);
  testMethod(testPetCompactThinPlateSplinesC2);
  testMethod(testPetCompactPolynomialC0);
  testMethod(testPetCompactPolynomialC6);
}

void PetRadialBasisFctMappingTest:: testPetThinPlateSplines()
{
  preciceTrace ( "testPetThinPlateSplines" );
  PetThinPlateSplines fct;
  PetRadialBasisFctMapping<PetThinPlateSplines> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<PetThinPlateSplines> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<PetThinPlateSplines> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<PetThinPlateSplines> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetMultiquadrics()
{
  preciceTrace ( "testPetMultiquadrics" );
  PetMultiquadrics fct(1e-3);
  PetRadialBasisFctMapping<PetMultiquadrics> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<PetMultiquadrics> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<PetMultiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<PetMultiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetInverseMultiquadrics()
{
  preciceTrace ( "testInverseMultiquadrics" );
  PetInverseMultiquadrics fct(1e-3);
  PetRadialBasisFctMapping<PetInverseMultiquadrics> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<PetInverseMultiquadrics> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<PetInverseMultiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<PetInverseMultiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetVolumeSplines()
{
  preciceTrace ( "testVolumeSplines" );
  PetVolumeSplines fct;
  PetRadialBasisFctMapping<PetVolumeSplines> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<PetVolumeSplines> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<PetVolumeSplines> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<PetVolumeSplines> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetGaussian()
{
  preciceTrace ( "testGaussian" );
  PetGaussian fct(1.0);
  PetRadialBasisFctMapping<PetGaussian> consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<PetGaussian> consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<PetGaussian> conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<PetGaussian> conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetCompactThinPlateSplinesC2()
{
  preciceTrace ( "testCompactThinPlateSplinesC2" );
  double supportRadius = 1.2;
  PetCompactThinPlateSplinesC2 fct(supportRadius);
  typedef PetRadialBasisFctMapping<PetCompactThinPlateSplinesC2> Mapping;
  Mapping consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetCompactPolynomialC0()
{
  preciceTrace ( "testCompactPolynomialC0" );
  double supportRadius = 1.2;
  PetCompactPolynomialC0 fct(supportRadius);
  typedef PetRadialBasisFctMapping<PetCompactPolynomialC0> Mapping;
  Mapping consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetCompactPolynomialC6()
{
  preciceTrace ( "testCompactPolynomialC6" );
  double supportRadius = 1.2;
  PetCompactPolynomialC6 fct(supportRadius);
  typedef PetRadialBasisFctMapping<PetCompactPolynomialC6> Mapping;
  Mapping consistentMap2D(Mapping::CONSISTENT, fct);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, fct);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, fct);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, fct);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: perform2DTestConsistentMapping
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
  validateNumericalEqualsWithEps ( value, 1.0, tolerance );

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
  validateNumericalEqualsWithEps ( value, 1.0, tolerance );

  vertex.setCoords ( Vector2D(1.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEqualsWithEps ( value, 2.0, tolerance );

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
  validateNumericalEqualsWithEps ( value, 2.0, tolerance );

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

void PetRadialBasisFctMappingTest:: perform2DTestConservativeMapping
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
  validate ( equals(values, tarch::la::Vector<4,double>(0.5, 0.5, 1.0, 1.0), tolerance) );

  vertex0.setCoords ( Vector2D(0.0, 0.5) );
  vertex1.setCoords ( Vector2D(1.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(values, tarch::la::Vector<4,double>(0.5, 1.0, 1.0, 0.5), tolerance) );

  vertex0.setCoords ( Vector2D(0.0, 1.0) );
  vertex1.setCoords ( Vector2D(1.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(values, tarch::la::Vector<4,double>(0.0, 2.0, 0.0, 1.0), tolerance) );
  
  vertex0.setCoords ( Vector2D(0.0, 0.0) );
  vertex1.setCoords ( Vector2D(1.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(values, tarch::la::Vector<4,double>(1.0, 0.0, 2.0, 0.0), tolerance) );

  vertex0.setCoords ( Vector2D(0.4, 0.5) );
  vertex1.setCoords ( Vector2D(0.6, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( sum(values), 3.0 );
}

void PetRadialBasisFctMappingTest:: perform3DTestConsistentMapping
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
  validateNumericalEqualsWithEps(value, 1.0, tolerance);

  vertex.setCoords(Vector3D(0.0, 0.5, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEqualsWithEps(value, 1.0, tolerance);

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
  validateNumericalEqualsWithEps(value, 1.0, tolerance);

  vertex.setCoords(Vector3D(1.0, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEqualsWithEps(value, 1.0, tolerance);

  vertex.setCoords(Vector3D(0.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEqualsWithEps(value, 2.0, tolerance);

  vertex.setCoords(Vector3D(1.0, 0.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEqualsWithEps(value, 2.0, tolerance);

  vertex.setCoords(Vector3D(1.0, 1.0, 1.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEqualsWithEps(value, 2.0, tolerance);

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
  validateNumericalEqualsWithEps(value, 1.5, tolerance);

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
  validateNumericalEqualsWithEps(value, 1.5, tolerance);

  vertex.setCoords(Vector3D(0.5, 0.5, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEquals(value, 1.5);
}

void PetRadialBasisFctMappingTest:: perform3DTestConservativeMapping
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
  validateWithParams1(tarch::la::equals(sum(values), expectedSum, tolerance), values);

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
