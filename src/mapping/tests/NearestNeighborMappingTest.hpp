// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MAPPING_NEARESTNEIGHBORMAPPINGTEST_HPP_
#define PRECICE_MAPPING_NEARESTNEIGHBORMAPPINGTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace mapping {
namespace tests {

/**
 * @brief Provides tests for class NearestNeighborMapping.
 */
class NearestNeighborMappingTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  NearestNeighborMappingTest();

  /**
   * @brief Destructor, empty.
   */
  virtual ~NearestNeighborMappingTest() {}

  /**
   * @brief Prepares run of tests, empty.
   */
  virtual void setUp() {}

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  void testConsistentNonIncremental();

  void testConservativeNonIncremental();
};

}}} // namespace precice, mapping, tests

#endif /* PRECICE_MAPPING_NEARESTNEIGHBORMAPPINGTEST_HPP_ */
