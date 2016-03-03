#ifndef PRECICE_NO_PETSC

#include "PetRadialBasisFctMappingTest.hpp"
#include <Eigen/Core>
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
  tolerance(1e-6)
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
      testMethod(testDistributedConsistent2DV1);
      testMethod(testDistributedConsistent2DV2);
      testMethod(testDistributedConservative2DV1);
      testMethod(testDistributedConservative2DV2);
      // testMethod(testDistributedConservative2DV3);
      // testMethod(testDistributedConservative2DV4);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
}

/// Test with a homogenous distribution of mesh amoung ranks
void PetRadialBasisFctMappingTest::testDistributedConsistent2DV1()
{
  preciceTrace("testDistributedConsistent2DV1");
  assertion(utils::Parallel::getCommunicatorSize() == 4);
  Gaussian fct(5.0);
  PetRadialBasisFctMapping<Gaussian> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);
  
  testDistributed(mapping,
                  { // Consistent mapping: The inMesh is communicated
                    {-1, 0, {0, 0}, {1}},
                    {-1, 0, {0, 1}, {2}},
                    {-1, 1, {1, 0}, {3}},
                    {-1, 1, {1, 1}, {4}},
                    {-1, 2, {2, 0}, {5}},
                    {-1, 2, {2, 1}, {6}},
                    {-1, 3, {3, 0}, {7}},
                    {-1, 3, {3, 1}, {8}}
                  },
                  { // The outMesh is local, every rank is used
                    {0, -1, {0, 0}, {0}},
                    {0, -1, {0, 1}, {0}},
                    {1, -1, {1, 0}, {0}},
                    {1, -1, {1, 1}, {0}},
                    {2, -1, {2, 0}, {0}},
                    {2, -1, {2, 1}, {0}},
                    {3, -1, {3, 0}, {0}},
                    {3, -1, {3, 1}, {0}}
                  },
                  { // Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                    { 0, {1} },
                    { 0, {2} },
                    { 1, {3} },
                    { 1, {4} },
                    { 2, {5} },
                    { 2, {6} },
                    { 3, {7} },
                    { 3, {8} }
                  }
    );
}

/// Using a more heterogenous distributon of vertices and owner
void PetRadialBasisFctMappingTest::testDistributedConsistent2DV2()
{
  preciceTrace("testDistributedConsistent2DV2");
  assertion(utils::Parallel::getCommunicatorSize() == 4);
  Gaussian fct(5.0);
  PetRadialBasisFctMapping<Gaussian> mapping(Mapping::CONSISTENT, 2, fct, false, false, false);
  
  testDistributed(mapping,
                  { // Consistent mapping: The inMesh is communicated, rank 2 owns no vertices
                    {-1, 0, {0, 0}, {1}},
                    {-1, 0, {0, 1}, {2}},
                    {-1, 1, {1, 0}, {3}},
                    {-1, 1, {1, 1}, {4}},
                    {-1, 1, {2, 0}, {5}},
                    {-1, 3, {2, 1}, {6}},
                    {-1, 3, {3, 0}, {7}},
                    {-1, 3, {3, 1}, {8}}
                  },
                  { // The outMesh is local, rank 1 is empty
                    {0, -1, {0, 0}, {0}},
                    {0, -1, {0, 1}, {0}},
                    {0, -1, {1, 0}, {0}},
                    {2, -1, {1, 1}, {0}},
                    {2, -1, {2, 0}, {0}},
                    {2, -1, {2, 1}, {0}},
                    {3, -1, {3, 0}, {0}},
                    {3, -1, {3, 1}, {0}}
                  },
                  { // Tests for {0, 1, 2} on the first rank,
                    // second rank (consistent with the outMesh) is empty, ...
                    { 0, {1} },
                    { 0, {2} },
                    { 0, {3} },
                    { 2, {4} },
                    { 2, {5} },
                    { 2, {6} },
                    { 3, {7} },
                    { 3, {8} }
                  }
    );
}

/// Test with a homogenous distribution of mesh amoung ranks
void PetRadialBasisFctMappingTest::testDistributedConservative2DV1()
{
  preciceTrace("testDistributedConservative2DV1");
  assertion(utils::Parallel::getCommunicatorSize() == 4);
  Gaussian fct(5.0);
  PetRadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);
  
  testDistributed(mapping,
                  { // Conservative mapping: The inMesh is local
                    {0, -1, {0, 0}, {1}},
                    {0, -1, {0, 1}, {2}},
                    {1, -1, {1, 0}, {3}},
                    {1, -1, {1, 1}, {4}},
                    {2, -1, {2, 0}, {5}},
                    {2, -1, {2, 1}, {6}},
                    {3, -1, {3, 0}, {7}},
                    {3, -1, {3, 1}, {8}}
                  },
                  { // The outMesh is distributed
                    {-1, 0, {0, 0}, {0}},
                    {-1, 0, {0, 1}, {0}},
                    {-1, 1, {1, 0}, {0}},
                    {-1, 1, {1, 1}, {0}},
                    {-1, 2, {2, 0}, {0}},
                    {-1, 2, {2, 1}, {0}},
                    {-1, 3, {3, 0}, {0}},
                    {-1, 3, {3, 1}, {0}}
                  },
                  { // Tests for {0, 1, 0, 0, 0, 0, 0, 0} on the first rank,
                    // {0, 0, 2, 3, 0, 0, 0, 0} on the second, ...
                    {0, {1}}, {0, {2}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}},
                    {1, {0}}, {1, {0}}, {1, {3}}, {1, {4}}, {1, {0}}, {1, {0}}, {1, {0}}, {1, {0}},
                    {2, {0}}, {2, {0}}, {2, {0}}, {2, {0}}, {2, {5}}, {2, {6}}, {2, {0}}, {2, {0}},
                    {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {7}}, {3, {8}}
                  },
                  utils::Parallel::getProcessRank()*2
    );
}

/// Using a more heterogenous distribution of vertices and owner
void PetRadialBasisFctMappingTest::testDistributedConservative2DV2()
{
  preciceTrace("testDistributedConservative2DV2");
  assertion(utils::Parallel::getCommunicatorSize() == 4);
  Gaussian fct(5.0);
  PetRadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 0, 4, 6};
  
  testDistributed(mapping,
                  { // Conservative mapping: The inMesh is local but rank 0 has no vertices
                    {1, -1, {0, 0}, {1}},
                    {1, -1, {0, 1}, {2}},
                    {1, -1, {1, 0}, {3}},
                    {1, -1, {1, 1}, {4}},
                    {2, -1, {2, 0}, {5}},
                    {2, -1, {2, 1}, {6}},
                    {3, -1, {3, 0}, {7}},
                    {3, -1, {3, 1}, {8}}
                  },
                  { // The outMesh is distributed, rank 0 owns no vertex
                    {-1, 1, {0, 0}, {0}},
                    {-1, 1, {0, 1}, {0}},
                    {-1, 1, {1, 0}, {0}},
                    {-1, 1, {1, 1}, {0}},
                    {-1, 2, {2, 0}, {0}},
                    {-1, 2, {2, 1}, {0}},
                    {-1, 3, {3, 0}, {0}},
                    {-1, 3, {3, 1}, {0}}
                  },
                  { // Tests for {0, 0, 0, 0, 0, 0, 0, 0} on the first rank,
                    // {1, 2, 2, 3, 0, 0, 0, 0} on the second, ...
                    {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}},
                    {1, {1}}, {1, {2}}, {1, {3}}, {1, {4}}, {1, {0}}, {1, {0}}, {1, {0}}, {1, {0}},
                    {2, {0}}, {2, {0}}, {2, {0}}, {2, {0}}, {2, {5}}, {2, {6}}, {2, {0}}, {2, {0}},
                    {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {7}}, {3, {8}}
                  },
                  globalIndexOffsets[utils::Parallel::getProcessRank()]
    );
}

/// Using meshes of different sizes, inMesh is smaller then outMesh
void PetRadialBasisFctMappingTest::testDistributedConservative2DV3()
{
  preciceTrace("testDistributedConservative2DV3");
  assertion(utils::Parallel::getCommunicatorSize() == 4);
  Gaussian fct(2.0);
  PetRadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 0, 3, 5};
  
  testDistributed(mapping,
                  { // Conservative mapping: The inMesh is local but rank 0 has no vertices
                    {1, -1, {0, 0}, {1}},
                    {1, -1, {1, 0}, {3}},
                    {1, -1, {1, 1}, {4}},
                    {2, -1, {2, 0}, {5}},
                    {2, -1, {2, 1}, {6}},
                    {3, -1, {3, 0}, {7}},
                    {3, -1, {3, 1}, {8}}
                  },
                  { // The outMesh is distributed, rank 0 owns no vertex
                    {-1, 1, {0, 0}, {0}},
                    {-1, 1, {0, 1}, {0}},
                    {-1, 1, {1, 0}, {0}},
                    {-1, 1, {1, 1}, {0}},
                    {-1, 2, {2, 0}, {0}},
                    {-1, 2, {2, 1}, {0}},
                    {-1, 3, {3, 0}, {0}},
                    {-1, 3, {3, 1}, {0}}
                  },
                  { // Tests for {0, 0, 0, 0, 0, 0, 0, 0} on the first rank,
                    // {1, 2, 2, 3, 0, 0, 0, 0} on the second, ...
                    {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}},
                    {1, {1}}, {1, {0}}, {1, {3}}, {1, {4}}, {1, {0}}, {1, {0}}, {1, {0}}, {1, {0}},
                    {2, {0}}, {2, {0}}, {2, {0}}, {2, {0}}, {2, {5}}, {2, {0}}, {2, {0}}, {2, {0}},
                    {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {7}}, {3, {8}}
                  },
                  globalIndexOffsets[utils::Parallel::getProcessRank()]
    );
}

/// Using meshes of different sizes, outMesh is smaller then inMesh
void PetRadialBasisFctMappingTest::testDistributedConservative2DV4()
{
  preciceTrace("testDistributedConservative2DV4");
  assertion(utils::Parallel::getCommunicatorSize() == 4);
  Gaussian fct(4.0);
  PetRadialBasisFctMapping<Gaussian> mapping(Mapping::CONSERVATIVE, 2, fct, false, false, false);

  std::vector<int> globalIndexOffsets = {0, 2, 4, 6};
  
  testDistributed(mapping,
                  { // Conservative mapping: The inMesh is local
                    {0, -1, {0, 0}, {1}},
                    {0, -1, {0, 1}, {2}},
                    {1, -1, {1, 0}, {3}},
                    {1, -1, {1, 1}, {4}},
                    {2, -1, {2, 0}, {5}},
                    {2, -1, {2, 1}, {6}},
                    {3, -1, {3, 0}, {7}},
                    {3, -1, {3, 1}, {8}}
                  },
                  { // The outMesh is distributed, rank has no vertex at all
                    {-1, 1, {0, 1}, {0}},
                    {-1, 1, {1, 0}, {0}},
                    {-1, 1, {1, 1}, {0}},
                    {-1, 2, {2, 0}, {0}},
                    {-1, 2, {2, 1}, {0}},
                    {-1, 3, {3, 0}, {0}},
                    {-1, 3, {3, 1}, {0}}
                  },
                  { // Tests for {0, 0, 0, 0, 0, 0, 0, 0} on the first rank,
                    // {2, 3, 4, 3, 0, 0, 0, 0} on the second, ...
                    {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}}, {0, {0}},
                    {1, {2.42855}}, {1, {3.61905}}, {1, {4.14286}}, {1, {0}}, {1, {0}}, {1, {0}}, {1, {0}},
                    {2, {0}}, {2, {0}}, {2, {0}}, {2, {5.33335942629867876263}}, {2, {5.85714}}, {2, {0}}, {2, {0}},
                    {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {0}}, {3, {7.04763872499617693990}}, {3, {7.57143}}
                  },
                  globalIndexOffsets[utils::Parallel::getProcessRank()]
    );
}
// Python results:
// 2.42857  3.61905  4.14286  5.33333  5.85714  7.04762  7.57143


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

void PetRadialBasisFctMappingTest:: perform2DTestConsistentMapping(Mapping& mapping)
{
  preciceTrace ( "perform2DTestConsistentMapping()" );
  int dimensions = 2;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  inMesh->createVertex ( Vector2d(0.0, 0.0) );
  inMesh->createVertex ( Vector2d(1.0, 0.0) );
  inMesh->createVertex ( Vector2d(1.0, 1.0) );
  inMesh->createVertex ( Vector2d(0.0, 1.0) );
  inMesh->allocateDataValues ();
  addGlobalIndex(inMesh);
  
  auto& values = inData->values();
  values << 1.0, 2.0, 2.0, 1.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();
  mesh::Vertex& vertex = outMesh->createVertex ( Vector2d(0, 0) );
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes ( inMesh, outMesh );
  validateEquals ( mapping.hasComputedMapping(), false );

  vertex.setCoords ( Vector2d(0.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  double value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEqualsWithEps ( value, 1.0, tolerance );

  vertex.setCoords ( Vector2d(0.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.0 );

  vertex.setCoords ( Vector2d(0.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEqualsWithEps ( value, 1.0, tolerance );

  vertex.setCoords ( Vector2d(1.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEqualsWithEps ( value, 2.0, tolerance );

  vertex.setCoords ( Vector2d(1.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 2.0 );

  vertex.setCoords ( Vector2d(1.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEqualsWithEps ( value, 2.0, tolerance );

  vertex.setCoords ( Vector2d(0.5, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.5 );

  vertex.setCoords ( Vector2d(0.5, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.5 );

  vertex.setCoords ( Vector2d(0.5, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.5 );
}


void PetRadialBasisFctMappingTest:: perform2DTestConservativeMapping(Mapping& mapping)
{
  preciceTrace ( "perform2DTestConservativeMapping()" );
  int dimensions = 2;
  using Eigen::Vector2d;

  // Create mesh to map from
  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", dimensions, false) );
  mesh::PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  mesh::Vertex& vertex0 = inMesh->createVertex ( Vector2d(0,0) );
  mesh::Vertex& vertex1 = inMesh->createVertex ( Vector2d(0,0) );
  inMesh->allocateDataValues ();
  inData->values() << 1.0, 2.0;
  addGlobalIndex(inMesh);

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID ();
  outMesh->createVertex ( Vector2d(0.0, 0.0) );
  outMesh->createVertex ( Vector2d(1.0, 0.0) );
  outMesh->createVertex ( Vector2d(1.0, 1.0) );
  outMesh->createVertex ( Vector2d(0.0, 1.0) );
  outMesh->allocateDataValues ();
  addGlobalIndex(outMesh);

  auto& values = outData->values();

  mapping.setMeshes ( inMesh, outMesh );
  validateEquals ( mapping.hasComputedMapping(), false );

  vertex0.setCoords ( Vector2d(0.5, 0.0) );
  vertex1.setCoords ( Vector2d(0.5, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(utils::DynVector(values), tarch::la::Vector<4,double>(0.5, 0.5, 1.0, 1.0), tolerance) );

  vertex0.setCoords ( Vector2d(0.0, 0.5) );
  vertex1.setCoords ( Vector2d(1.0, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(utils::DynVector(values), tarch::la::Vector<4,double>(0.5, 1.0, 1.0, 0.5), tolerance) );

  vertex0.setCoords ( Vector2d(0.0, 1.0) );
  vertex1.setCoords ( Vector2d(1.0, 0.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(utils::DynVector(values), tarch::la::Vector<4,double>(0.0, 2.0, 0.0, 1.0), tolerance) );
  
  vertex0.setCoords ( Vector2d(0.0, 0.0) );
  vertex1.setCoords ( Vector2d(1.0, 1.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validate ( equals(utils::DynVector(values), tarch::la::Vector<4,double>(1.0, 0.0, 2.0, 0.0), tolerance) );

  vertex0.setCoords ( Vector2d(0.4, 0.5) );
  vertex1.setCoords ( Vector2d(0.6, 0.5) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( values.sum(), 3.0 );
}

void PetRadialBasisFctMappingTest:: perform3DTestConsistentMapping(Mapping& mapping)
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
  
  auto& values = inData->values();
  values << 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0;

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

void PetRadialBasisFctMappingTest:: perform3DTestConservativeMapping(Mapping& mapping)
{
  preciceTrace ( "perform3DTestConservativeMapping()" );
  using Eigen::Vector3d;
  int dimensions = 3;

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, false));
  mesh::PtrData inData = inMesh->createData("InData", 1);
  int inDataID = inData->getID();
  mesh::Vertex& vertex0 = inMesh->createVertex(Vector3d(0,0,0));
  mesh::Vertex& vertex1 = inMesh->createVertex(Vector3d(0,0,0));
  inMesh->allocateDataValues();
  inData->values() << 1.0, 2.0;
  addGlobalIndex(inMesh);
  
  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, false));
  mesh::PtrData outData = outMesh->createData("OutData", 1);
  int outDataID = outData->getID();
  outMesh->createVertex(Vector3d(0.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 0.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 0.0));
  outMesh->createVertex(Vector3d(0.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 0.0, 1.0));
  outMesh->createVertex(Vector3d(1.0, 1.0, 1.0));
  outMesh->createVertex(Vector3d(0.0, 1.0, 1.0));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  
  auto& values = outData->values();
  double expectedSum = inData->values().sum();

  mapping.setMeshes(inMesh, outMesh);
  validateEquals(mapping.hasComputedMapping(), false);

  vertex0.setCoords(Vector3d(0.5, 0.0, 0.0));
  vertex1.setCoords(Vector3d(0.5, 1.0, 0.0));
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  validateEquals(mapping.hasComputedMapping(), true);
  validateWithParams1(tarch::la::equals(values.sum(), expectedSum, tolerance), values);
}

void PetRadialBasisFctMappingTest:: testDeadAxis2D()
{
  preciceTrace ( "testDeadAxis2D()" );
  using Eigen::Vector2d;
  int dimensions = 2;
  
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
  inMesh->createVertex ( Vector2d(0.0, 1.0) );
  inMesh->createVertex ( Vector2d(1.0, 1.0) );
  inMesh->createVertex ( Vector2d(2.0, 1.0) );
  inMesh->createVertex ( Vector2d(3.0, 1.0) );
  inMesh->allocateDataValues ();
  addGlobalIndex(inMesh);
  
  auto& values = inData->values();
  values << 1.0, 2.0, 2.0, 1.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();
  mesh::Vertex& vertex = outMesh->createVertex ( Vector2d(0,0) );
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes ( inMesh, outMesh );
  validateEquals ( mapping.hasComputedMapping(), false );

  vertex.setCoords ( Vector2d(0.0, 3.0) );
  mapping.computeMapping ();
  mapping.map ( inDataID, outDataID );
  double value = outData->values()[0];
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( value, 1.0 );
}

void PetRadialBasisFctMappingTest:: testDeadAxis3D()
{
  preciceTrace ( "testDeadAxis3D()" );
  using Eigen::Vector3d;
  int dimensions = 3;

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
  inMesh->createVertex ( Vector3d(0.0, 3.0, 0.0) );
  inMesh->createVertex ( Vector3d(1.0, 3.0, 0.0) );
  inMesh->createVertex ( Vector3d(0.0, 3.0, 1.0) );
  inMesh->createVertex ( Vector3d(1.0, 3.0, 1.0) );
  inMesh->allocateDataValues ();
  addGlobalIndex(inMesh);
  
  auto& values = inData->values();
  values << 1.0, 2.0, 3.0, 4.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh ( new mesh::Mesh("OutMesh", dimensions, false) );
  mesh::PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();
  outMesh->createVertex ( Vector3d(0.0, 2.9, 0.0) );
  outMesh->createVertex ( Vector3d(0.8, 2.9, 0.1) );
  outMesh->createVertex ( Vector3d(0.1, 2.9, 0.9) );
  outMesh->createVertex ( Vector3d(1.1, 2.9, 1.1) );
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

/*
MeshSpecification format:
{ {rank, owner rank, {x, y, z}, {v}}, ... }

also see struct VertexSpecification at the header.

- -1 on rank means all ranks
- -1 on owner rank means no rank
- x, y, z is position of vertex, z is optional, 2D mesh will be created then
- v is the value of the respective vertex. Only 1D supported at this time.

ReferenceSpecification format:
{ {rank, {v}, ... }
- -1 on rank means all ranks
- v is the expected value of n-th vertex on that particular rank
*/

void PetRadialBasisFctMappingTest::getDistributedMesh(MeshSpecification const & vertices,
                                                      mesh::PtrMesh& mesh,
                                                      mesh::PtrData& data,
                                                      int globalIndexOffset)
{
  preciceTrace("getDistributedMesh");
  using Par = utils::Parallel;
  Eigen::VectorXd d;

  int i = 0;
  for (auto& vertex : vertices) {
    if (vertex.rank == Par::getProcessRank() or vertex.rank == -1) {
      if (vertex.position.size() == 3)  // 3-dimensional
        mesh->createVertex(Eigen::Vector3d(vertex.position.data()));
      else if (vertex.position.size() == 2) // 2-dimensional
        mesh->createVertex(Eigen::Vector2d(vertex.position.data()));
      
      if (vertex.owner == Par::getProcessRank())
        mesh->vertices().back().setOwner(true);
      else
        mesh->vertices().back().setOwner(false);

      d.conservativeResize(i+1);
      d[i] = vertex.value[0]; // only 1-d value here for now
      i++;
    }

  }
  addGlobalIndex(mesh, globalIndexOffset );
  mesh->allocateDataValues();
  data->values() = d;
}


void PetRadialBasisFctMappingTest::testDistributed(Mapping& mapping,
                                                   MeshSpecification inMeshSpec,
                                                   MeshSpecification outMeshSpec,
                                                   ReferenceSpecification referenceSpec,
                                                   int inGlobalIndexOffset)                                        {
  preciceTrace("testDistributed");
  using Par = utils::Parallel;
  assertion(Par::getCommunicatorSize() == 4);
  int meshDimension = inMeshSpec[0].position.size();
  int valueDimension = inMeshSpec[0].value.size();

  mesh::PtrMesh inMesh ( new mesh::Mesh("InMesh", meshDimension, false) );
  mesh::PtrData inData = inMesh->createData("InData", valueDimension);
  int inDataID = inData->getID();

  getDistributedMesh(inMeshSpec, inMesh, inData, inGlobalIndexOffset);

  mesh::PtrMesh outMesh ( new mesh::Mesh("outMesh", meshDimension, false) );
  mesh::PtrData outData = outMesh->createData( "OutData", valueDimension );
  int outDataID = outData->getID();

  getDistributedMesh(outMeshSpec, outMesh, outData);
  
  mapping.setMeshes(inMesh, outMesh);
  validateEquals(mapping.hasComputedMapping(), false);

  mapping.computeMapping();
  validateEquals(mapping.hasComputedMapping(), true);
  mapping.map(inDataID, outDataID);

  int index = 0;
  for (auto& referenceVertex : referenceSpec) {
    if (referenceVertex.first == Par::getProcessRank() or referenceVertex.first == -1) {
      for (auto& point : referenceVertex.second) {
        // only 1-d here for now
        validateNumericalEqualsWithEps(outData->values()[index], point, tolerance);
      }
      ++index;
    }
  }
  validateEquals(outData->values().size(), index);
}



}}} // namespace precice, mapping, tests

#endif
