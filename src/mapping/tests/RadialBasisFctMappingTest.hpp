// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MAPPING_TESTS_RADIALBASISFCTMAPPINGTEST_HPP_
#define PRECICE_MAPPING_TESTS_RADIALBASISFCTMAPPINGTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

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
 * @brief Provides tests for class RadialBasisFctMapping.
 */
class RadialBasisFctMappingTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  RadialBasisFctMappingTest();

  /**
   * @brief Destructor, empty.
   */
  virtual ~RadialBasisFctMappingTest() {}

  /**
   * @brief Prepares test run, empty.
   */
  virtual void setUp() {}

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  // @brief Logging device.
  static logging::Logger _log;

  void testThinPlateSplines();

  void testMultiquadrics();

  void testInverseMultiquadrics();

  void testVolumeSplines();

  void testGaussian();

  void testCompactThinPlateSplinesC2();

  void testCompactPolynomialC0();

  void testCompactPolynomialC6();

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

  /**
   * @brief
   */
  void testDeadAxis2D ();

  /**
   * @brief
   */
  void testDeadAxis3D ();
};

}}} // namespace precice, mapping, tests

#endif /* PRECICE_MAPPING_TESTS_RADIALBASISFCTMAPPINGTEST_HPP_ */
