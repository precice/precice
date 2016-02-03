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
  
  using Par = utils::Parallel;
  if (Par::getCommunicatorSize() > 3) {
    MPI_Comm comm = Par::getRestrictedCommunicator( {0, 1, 2, 3} );
    if (Par::getProcessRank() <= 3){
      Par::setGlobalCommunicator(comm); // hier auch noch PETSC_COMM_WORLD neu setzen?
      testMethod(testDistributedConsistent2D);
      testMethod(testDistributedConservative2D);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
}

/*
 * Creates a mesh that lies on different ranks. The list of vertices is formatted like that:
 * { {rank, owner rank, value, x, y, z}, {...}, ... }
 * If rank == -1, the vertex lives on all ranks
 * If the z value is ommited, a 2-dimensional mesh is generated.
 */
void PetRadialBasisFctMappingTest::getDistributedMesh(const std::vector<std::vector<int>> vertices,
                                                    mesh::PtrMesh& mesh,
                                                    mesh::PtrData& data)
{
  preciceTrace("getDistributedMesh");
  using Par = utils::Parallel;
  utils::DynVector d;
  
  for (auto& vertex : vertices) {
    if (Par::getProcessRank() == vertex[0] or vertex[0] == -1) {
      if (vertex.size() == 6)  // 3-dimensional
        mesh->createVertex(utils::Vector3D(vertex[3], vertex[4], vertex[5]));
      else if (vertex.size() == 5) // 2-dimensional
        mesh->createVertex(utils::Vector2D(vertex[3], vertex[4]));
      
      if (vertex[1] >= 0 and Par::getProcessRank() == vertex[1])
        mesh->vertices().back().setOwner(true);
      else
        mesh->vertices().back().setOwner(false);
      
      d.append(vertex[2]);
    }

  }
  addGlobalIndex(mesh);
  mesh->allocateDataValues();
  preciceDebug("Data values to be assigned = " << d);
  data->values() = d;
}


void PetRadialBasisFctMappingTest::testDistributedConsistent2D()
{
  preciceTrace("testDistributedConsistent2D");
  using Par = utils::Parallel;
  using utils::Vector2D;
  assertion(Par::getCommunicatorSize() == 4);
  int dimensions = 2;

  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", 2, false) );
  mesh::PtrData inData = inMesh->createData( "InData", 1 );
  int inDataID = inData->getID();

  // Consistent mapping: The inMesh is communicated, the outMesh is local
  getDistributedMesh(
    {
      {-1, 0, 0, 0, 0},
      {-1, 0, 1, 0, 1},
      {-1, 1, 2, 1, 0},
      {-1, 1, 3, 1, 1},
      {-1, 2, 4, 2, 0},
      {-1, 2, 5, 2, 1},
      {-1, 3, 6, 3, 0},
      {-1, 3, 7, 3, 1}
    }, inMesh, inData);

  mesh::PtrMesh outMesh ( new mesh::Mesh("outMesh", 2, false) );
  mesh::PtrData outData = outMesh->createData( "OutData", 1 );
  int outDataID = outData->getID();

  getDistributedMesh(
    {
      {0, -1, 0, 0, 0},
      {0, -1, 0, 0, 1},
      {1, -1, 0, 1, 0},
      {1, -1, 0, 1, 1},
      {2, -1, 0, 2, 0},
      {2, -1, 0, 2, 1},
      {3, -1, 0, 3, 0},
      {3, -1, 0, 3, 1}
    }, outMesh, outData);
  
  Gaussian fct(2.0);
  PetRadialBasisFctMapping<Gaussian> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);
  mapping.setMeshes(inMesh, outMesh);
  validateEquals(mapping.hasComputedMapping(), false);

  mapping.computeMapping();
  validateEquals(mapping.hasComputedMapping(), true);
  mapping.map(inDataID, outDataID);

  // Tests for {0, 1} on the first rank, {1, 2} on the second, ...
  for (size_t i=0; i < outData->values().size(); i++) {
    preciceDebug("outData->values()[" << i <<  "] = " << outData->values()[i]);
    validateNumericalEqualsWithEps ( outData->values()[i], Par::getProcessRank()*2 + i, tolerance );
  }
 
  // inMesh->computeState();
  // inMesh->computeDistribution();
  // computeState && computeDistribution auf mesh aufrufen ODER
  // addGlobalID mit Offset aufrufen je nach Größe der Basisfunktion vertices zweimal haben, aber nur einmal mit owner(true)
  // Auch Testcase mit einem Prozessor leer lassen

  // for (mesh::Vertex& v : inMesh->vertices()) {
    // std::cout << "Vertex globalIndex = " << v.getGlobalIndex() << std::endl;
  // }
}


void PetRadialBasisFctMappingTest::testDistributedConservative2D()
{
  preciceTrace("testDistributedConservative2D");
  using Par = utils::Parallel;
  using utils::Vector2D;
  assertion(Par::getCommunicatorSize() == 4);
  int dimensions = 2;

  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", 2, false) );
  mesh::PtrData inData = inMesh->createData( "InData", 1 );
  int inDataID = inData->getID();

  // Conservative mapping: The inMesh is local, the outMesh is communicated
  getDistributedMesh(
    {
      {0, -1, 0, 0, 0},
      {0, -1, 1, 0, 1},
      {1, -1, 2, 1, 0},
      {1, -1, 3, 1, 1},
      {2, -1, 4, 2, 0},
      {2, -1, 5, 2, 1},
      {3, -1, 6, 3, 0},
      {3, -1, 7, 3, 1}
    }, inMesh, inData);

  addGlobalIndex(inMesh, Par::getProcessRank()*2);

  mesh::PtrMesh outMesh ( new mesh::Mesh("outMesh", 2, false) );
  mesh::PtrData outData = outMesh->createData( "OutData", 1 );
  int outDataID = outData->getID();

  getDistributedMesh(
    {
      {-1, 0, 0, 0, 0},
      {-1, 0, 0, 0, 1},
      {-1, 1, 0, 1, 0},
      {-1, 1, 0, 1, 1},
      {-1, 2, 0, 2, 0},
      {-1, 2, 0, 2, 1},
      {-1, 3, 0, 3, 0},
      {-1, 3, 0, 3, 1}
    }, outMesh, outData);
  
  // preciceDebug("inData = " << inData->values());
  // for (auto& v : inMesh->vertices()) {
  //   if (v.isOwner()) {
  //     preciceDebug("owned Vertex = " << v.getCoords() << " Value = " << inData->values()[v.getID()]); }
  //   else {
  //     preciceDebug("Not owned Vertex = " << v.getCoords() << " Value = " << inData->values()[v.getID()]); }
  // }

  Gaussian fct(2.0);
  PetRadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);
  mapping.setMeshes(inMesh, outMesh);
  validateEquals(mapping.hasComputedMapping(), false);
  mapping.computeMapping();
  validateEquals(mapping.hasComputedMapping(), true);
  mapping.map(inDataID, outDataID);
  

  preciceDebug("Final out values = " << outData->values());

  tarch::la::DynamicVector<double> reference(8, 0);
  reference[Par::getProcessRank()*2]     = Par::getProcessRank()*2;
  reference[Par::getProcessRank()*2 + 1] = Par::getProcessRank()*2 + 1;
  preciceDebug("Reference = " << reference);
  for (int i = 0; i < reference.size(); i++) {
    validateNumericalEqualsWithEps( outData->values()[i], reference[i], tolerance );
  }
  
  // validateNumericalVectorEquals(outData->values(), reference);
        
  // Tests for {0, 1} on the first rank, {1, 2} on the second, ...
  // for (size_t i=0; i < outData->values().size(); i++) {
    // preciceDebug("outData->values()[" << i <<  "] = " << outData->values()[i]);
    // validateNumericalEqualsWithEps ( outData->values()[i], Par::getProcessRank()*2 + i, tolerance );
  // }
 
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
