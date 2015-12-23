#ifndef PRECICE_NO_PETSC

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

tarch::logging::Log PetRadialBasisFctMappingTest::_log ("precice::mapping::tests::PetRadialBasisFctMappingTest");

PetRadialBasisFctMappingTest:: PetRadialBasisFctMappingTest()
  :
  TestCase ( "precice::mapping::tests::PetRadialBasisFctMappingTest" ),
  tolerance(1e-7)
{}

void PetRadialBasisFctMappingTest:: run()
{
  PRECICE_MASTER_ONLY {
    PETSC_COMM_WORLD = PETSC_COMM_SELF;
    testMethod(testPetThinPlateSplines);
    testMethod(testPetMultiquadrics);
    testMethod(testPetInverseMultiquadrics);
    testMethod(testPetVolumeSplines);
    testMethod(testPetGaussian);
    testMethod(testPetCompactThinPlateSplinesC2);
    testMethod(testPetCompactPolynomialC0);
    testMethod(testPetCompactPolynomialC6);
    testMethod(testDeadAxis2D);
    testMethod(testDeadAxis3D);
    PETSC_COMM_WORLD = MPI_COMM_WORLD;
  }
  
  // testMethod(testMPI);
  
  using Par = utils::Parallel;
  if (Par::getCommunicatorSize() > 3) {
    // const std::vector<int> ranksWanted = {0, 1, 2 , 3};
    // MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
    MPI_Comm comm = Par::getRestrictedCommunicator( {0, 1, 2, 3} );
    if (Par::getProcessRank() <= 3){
      Par::setGlobalCommunicator(comm); // hier auch noch PETSC_COMM_WORLD neu setzen?
      testMethod(testMPI);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
}

mesh::PtrMesh PetRadialBasisFctMappingTest::getDistributedMesh(const std::vector<int>& ownedVertices)
{
  const std::vector<double> Xcoords = {0, 1, 2, 3};
  const std::vector<double> Ycoords = {0, 1};
  
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", 2, false) );
  for (auto& x : Xcoords) {
    for (auto& y : Ycoords) {
      inMesh->createVertex(utils::Vector2D(x, y));
      inMesh->vertices().back().setOwner(false);
    }
  }
  for (auto& owned : ownedVertices) {
    inMesh->vertices()[owned].setOwner(true);
  }
  addGlobalIndex(inMesh);
  return inMesh;
}

void PetRadialBasisFctMappingTest::testMPI()
{
  preciceTrace("testMPI");
  using Par = utils::Parallel;
  // assertion(Par::getCommunicatorSize() == 4);
  using utils::Vector2D;
  int dimensions = 2;
  mesh::PtrMesh inMesh, outMesh;
  
  if (Par::getProcessRank() == 0) {
    inMesh = getDistributedMesh( {0, 1} );
    outMesh = getDistributedMesh( {0, 1} );
  }
  else if (Par::getProcessRank() == 1) {
    inMesh = getDistributedMesh( {2, 3} );
    outMesh = getDistributedMesh( {2, 3} );
  }
  else if (Par::getProcessRank() == 2) {
    inMesh = getDistributedMesh( {4, 5} );
    outMesh = getDistributedMesh( {4, 5} );
  }
  else if (Par::getProcessRank() == 3) {
    inMesh = getDistributedMesh( {6, 7} );
    outMesh = getDistributedMesh( {6, 7} );
  }
  mesh::PtrData inData = inMesh->createData( "InData", 1 );
  int inDataID = inData->getID();
  inMesh->allocateDataValues();

  mesh::PtrData outData = outMesh->createData( "OutData", 1 );
  int outDataID = outData->getID();
  outMesh->allocateDataValues();

  std::vector<int> vi = { 0, 1, 2, 3, 4, 5, 6, 7 };
  utils::DynVector& values = inData->values();
  values = vi;

  for (auto& v : inMesh->vertices()) {
    if (v.isOwner())
      preciceDebug("owned Vertex = " << v.getCoords());
  }

  Gaussian fct(2.0);
  PetRadialBasisFctMapping<Gaussian> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);
  mapping.setMeshes(inMesh, outMesh);
  validateEquals(mapping.hasComputedMapping(), false);

  mapping.computeMapping();
  // mapping.map(inDataID, outDataID);
  // validateEquals(mapping.hasComputedMapping(), true);
  
  
 
  // inMesh->computeState();
  // inMesh->computeDistribution();
  // computeState && computeDistribution auf mesh aufrufen ODER
  // addGlobalID mit Offset aufrufen je nach Größe der Basisfunktion vertices zweimal haben, aber nur einmal mit owner(true)
  // Auch Testcase mit einem Prozessor leer lassen

  // for (mesh::Vertex& v : inMesh->vertices()) {
    // std::cout << "Vertex globalIndex = " << v.getGlobalIndex() << std::endl;
  // }
 

}


void PetRadialBasisFctMappingTest:: testPetThinPlateSplines()
{
  preciceTrace ( "testPetThinPlateSplines" );
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  ThinPlateSplines fct;
  PetRadialBasisFctMapping<ThinPlateSplines> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<ThinPlateSplines> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<ThinPlateSplines> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<ThinPlateSplines> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetMultiquadrics()
{
  preciceTrace ( "testPetMultiquadrics" );
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  Multiquadrics fct(1e-3);
  PetRadialBasisFctMapping<Multiquadrics> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<Multiquadrics> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<Multiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<Multiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetInverseMultiquadrics()
{
  preciceTrace ( "testInverseMultiquadrics" );
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  InverseMultiquadrics fct(1e-3);
  PetRadialBasisFctMapping<InverseMultiquadrics> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<InverseMultiquadrics> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<InverseMultiquadrics> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<InverseMultiquadrics> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetVolumeSplines()
{
  preciceTrace ( "testVolumeSplines" );
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  VolumeSplines fct;
  PetRadialBasisFctMapping<VolumeSplines> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<VolumeSplines> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<VolumeSplines> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<VolumeSplines> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetGaussian()
{
  preciceTrace ( "testGaussian" );
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  Gaussian fct(1.0);
  PetRadialBasisFctMapping<Gaussian> consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  PetRadialBasisFctMapping<Gaussian> consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  PetRadialBasisFctMapping<Gaussian> conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  PetRadialBasisFctMapping<Gaussian> conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetCompactThinPlateSplinesC2()
{
  preciceTrace ( "testCompactThinPlateSplinesC2" );
  double supportRadius = 1.2;
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  CompactThinPlateSplinesC2 fct(supportRadius);
  typedef PetRadialBasisFctMapping<CompactThinPlateSplinesC2> Mapping;
  Mapping consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetCompactPolynomialC0()
{
  preciceTrace ( "testCompactPolynomialC0" );
  double supportRadius = 1.2;
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  CompactPolynomialC0 fct(supportRadius);
  typedef PetRadialBasisFctMapping<CompactPolynomialC0> Mapping;
  Mapping consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
  perform3DTestConservativeMapping(conservativeMap3D);
}

void PetRadialBasisFctMappingTest:: testPetCompactPolynomialC6()
{
  preciceTrace ( "testCompactPolynomialC6" );
  double supportRadius = 1.2;
  bool xDead = false;
  bool yDead = false;
  bool zDead = false;
  CompactPolynomialC6 fct(supportRadius);
  typedef PetRadialBasisFctMapping<CompactPolynomialC6> Mapping;
  Mapping consistentMap2D(Mapping::CONSISTENT, 2, fct, xDead, yDead, zDead);
  perform2DTestConsistentMapping(consistentMap2D);
  Mapping consistentMap3D(Mapping::CONSISTENT, 3, fct, xDead, yDead, zDead);
  perform3DTestConsistentMapping(consistentMap3D);
  Mapping conservativeMap2D(Mapping::CONSERVATIVE, 2, fct, xDead, yDead, zDead);
  perform2DTestConservativeMapping(conservativeMap2D);
  Mapping conservativeMap3D(Mapping::CONSERVATIVE, 3, fct, xDead, yDead, zDead);
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
  addGlobalIndex(inMesh);
  
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
  addGlobalIndex(outMesh);

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
  addGlobalIndex(inMesh);

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID ();
  outMesh->createVertex ( Vector2D(0.0, 0.0) );
  outMesh->createVertex ( Vector2D(1.0, 0.0) );
  outMesh->createVertex ( Vector2D(1.0, 1.0) );
  outMesh->createVertex ( Vector2D(0.0, 1.0) );
  outMesh->allocateDataValues ();
  addGlobalIndex(outMesh);

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
  addGlobalIndex(inMesh);
  
  utils::DynVector& values = inData->values();
  assignList(values) = 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false));
  mesh::PtrData outData = outMesh->createData("OutData", 1);
  int outDataID = outData->getID();
  mesh::Vertex& vertex = outMesh->createVertex(Vector3D(0.0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

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
  validateNumericalEqualsWithEps(value, 1.5, tolerance );

  vertex.setCoords(Vector3D(0.0, 1.0, 0.5));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  value = outData->values()[0];
  validateEquals(mapping.hasComputedMapping(), true);
  validateNumericalEqualsWithEps(value, 1.5, tolerance);

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
  addGlobalIndex(inMesh);
  
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
  addGlobalIndex(outMesh);
  
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

void PetRadialBasisFctMappingTest:: testDeadAxis2D
()
{
  preciceTrace ( "testDeadAxis2D()" );
  int dimensions = 2;
  using utils::Vector2D;

  bool xDead = false;
  bool yDead = true;
  bool zDead = false;

  ThinPlateSplines fct;
  PetRadialBasisFctMapping<ThinPlateSplines> mapping(Mapping::CONSISTENT, dimensions, fct,
                                                     xDead, yDead, zDead);

  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  inMesh->createVertex ( Vector2D(0.0, 1.0) );
  inMesh->createVertex ( Vector2D(1.0, 1.0) );
  inMesh->createVertex ( Vector2D(2.0, 1.0) );
  inMesh->createVertex ( Vector2D(3.0, 1.0) );
  inMesh->allocateDataValues ();
  addGlobalIndex(inMesh);
  
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
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes ( inMesh, outMesh );
  validateEquals ( mapping.hasComputedMapping(), false );

  vertex.setCoords ( Vector2D(0.0, 3.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  double value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.0 );
}

void PetRadialBasisFctMappingTest:: testDeadAxis3D()
{
  preciceTrace ( "testDeadAxis3D()" );
  int dimensions = 3;
  using utils::Vector3D;

  double supportRadius = 1.2;
  CompactPolynomialC6 fct(supportRadius);
  bool xDead = false;
  bool yDead = true;
  bool zDead = false;
  typedef PetRadialBasisFctMapping<CompactPolynomialC6> Mapping;
  Mapping mapping(Mapping::CONSISTENT, dimensions, fct, xDead, yDead, zDead);

  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  inMesh->createVertex ( Vector3D(0.0, 3.0, 0.0) );
  inMesh->createVertex ( Vector3D(1.0, 3.0, 0.0) );
  inMesh->createVertex ( Vector3D(0.0, 3.0, 1.0) );
  inMesh->createVertex ( Vector3D(1.0, 3.0, 1.0) );
  inMesh->allocateDataValues ();
  addGlobalIndex(inMesh);
  
  tarch::la::Vector<4,double> assignValues;
  assignList(assignValues) = 1.0, 2.0, 3.0, 4.0;
  utils::DynVector& values = inData->values();
  values = assignValues;

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();
  outMesh->createVertex ( Vector3D(0.0, 2.9, 0.0) );
  outMesh->createVertex ( Vector3D(0.8, 2.9, 0.1) );
  outMesh->createVertex ( Vector3D(0.1, 2.9, 0.9) );
  outMesh->createVertex ( Vector3D(1.1, 2.9, 1.1) );
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes ( inMesh, outMesh );
  validateEquals ( mapping.hasComputedMapping(), false );

  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );

  validateNumericalEquals ( outData->values()[0], 1.0 );
  validateNumericalEquals ( outData->values()[1], 2.0 );
  validateNumericalEquals ( outData->values()[2], 2.9 );
  validateNumericalEquals ( outData->values()[3], 4.3 );
}

void PetRadialBasisFctMappingTest::addGlobalIndex(mesh::PtrMesh &mesh, int offset)
{
  for (mesh::Vertex& v : mesh->vertices()) {
    v.setGlobalIndex(v.getID() + offset);
  }
}



}}} // namespace precice, mapping, tests

#endif
