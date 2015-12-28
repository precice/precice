// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MAPPING_TESTS_MAPPINGCONFIGURATIONTEST_HPP_
#define PRECICE_MAPPING_TESTS_MAPPINGCONFIGURATIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace mapping {
namespace tests {

/**
 * @brief Provides tests for class MappingConfiguration.
 */
class MappingConfigurationTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor, takes path to src directory of precice as argument.
   */
  MappingConfigurationTest ();

  /**
   * @brief Destructor, empty.
   */
  virtual ~MappingConfigurationTest () {};

  /**
   * @brief Retrieves path to test directory.
   */
  virtual void setUp ();

  /**
   * @brief Runs all tests.
   */
  virtual void run ();

private:

  // @brief Logging device.
  static logging::Logger _log;

  // @brief Path to src directory of precice.
  std::string _pathToTests;
};

}}} // namespace precice, mapping, tests

#endif /* PRECICE_MAPPING_TESTS_MAPPINGCONFIGURATIONTEST_HPP_ */
