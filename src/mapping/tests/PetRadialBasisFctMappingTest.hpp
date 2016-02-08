#ifndef PRECICE_NO_PETSC
#pragma once

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"
#include "mesh/Mesh.hpp"


namespace precice {
  namespace mapping {
    class Mapping;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mapping {
namespace tests {

/**
 * @brief Provides tests for class PetRadialBasisFctMapping.
 */
class PetRadialBasisFctMappingTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  PetRadialBasisFctMappingTest();

  /**
   * @brief Destructor, empty.
   */
  virtual ~PetRadialBasisFctMappingTest() {}

  /**
   * @brief Prepares test run, empty.
   */
  virtual void setUp() {}

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  /// Petsc uses iterative methods to solve the linear system.
  const double tolerance;

  // @brief Logging device.
  static tarch::logging::Log _log;

  void testDistributedConsistent2DV1();
  void testDistributedConsistent2DV2();

  void testDistributedConservative2D();

  void testPetThinPlateSplines();

  void testPetMultiquadrics();

  void testPetInverseMultiquadrics();

  void testPetVolumeSplines();

  void testPetGaussian();

  void testPetCompactThinPlateSplinesC2();

  void testPetCompactPolynomialC0();

  void testPetCompactPolynomialC6();

  void perform2DTestConsistentMapping ( Mapping& mapping );

  void perform2DTestConservativeMapping ( Mapping& mapping );

  void perform3DTestConsistentMapping ( Mapping& mapping );

  void perform3DTestConservativeMapping ( Mapping& mapping );

  void testDeadAxis2D();

  void testDeadAxis3D();

  struct VertexSpecification
  {
    int rank;
    int owner;
    std::vector<double> position;
    std::vector<double> value;
  };

  using MeshSpecification = std::vector<VertexSpecification>;

  /// Contains which values are expected on which rank: rank -> vector of data.
  using ReferenceSpecification = std::vector<std::pair<int,std::vector<double>>>;

  /// Helper function: Add the global index from vertex::getID
  void addGlobalIndex(mesh::PtrMesh &mesh, int offset = 0);

  /// Helper function: create a distributed mesh with values and owner set on some ranks.
  void getDistributedMesh(MeshSpecification const & vertices,
                          mesh::PtrMesh& mesh,
                          mesh::PtrData& data,
                          int globalIndexOffset = 0);
  

  /// Helper function for testing distributed meshes
  void testDistributed(Mapping& mapping,
                       MeshSpecification inMeshSpec,
                       MeshSpecification outMeshSpec,
                       ReferenceSpecification referenceSpec,
                       int inGlobalIndexOffset = 0);


};

}}} // namespace precice, mapping, tests


#endif
