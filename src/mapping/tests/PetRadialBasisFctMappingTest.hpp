#ifndef PRECICE_NO_PETSC
#pragma once

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

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

  /// @brief Petsc uses iterative methods to solve the linear system.
  const double tolerance;

  // @brief Logging device.
  static tarch::logging::Log _log;

  void testPetThinPlateSplines();

  void testPetMultiquadrics();

  void testPetInverseMultiquadrics();

  void testPetVolumeSplines();

  void testPetGaussian();

  void testPetCompactThinPlateSplinesC2();

  void testPetCompactPolynomialC0();

  void testPetCompactPolynomialC6();

  /**
   * @brief
   */
  void perform2DTestConsistentMapping ( Mapping& mapping );

  /**
   * @brief
   */
  void perform2DTestConservativeMapping ( Mapping& mapping );

  /**
   * @brief
   */
  void perform3DTestConsistentMapping ( Mapping& mapping );

  /**
   * @brief
   */
  void perform3DTestConservativeMapping ( Mapping& mapping );

  void testDeadAxis2D();

  void testDeadAxis3D();
};

}}} // namespace precice, mapping, tests


#endif
