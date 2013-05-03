// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_TESTS_HIERARCHICALAITKENPOSTPROCESSINGTEST_HPP_
#define PRECICE_CPLSCHEME_TESTS_HIERARCHICALAITKENPOSTPROCESSINGTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace cplscheme {
namespace tests {

/**
 * @brief Provides test cases for class HierarchicalAitkenPostProcessing.
 */
class HierarchicalAitkenPostProcessingTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  HierarchicalAitkenPostProcessingTest ();

  /**
   * @brief Destructor.
   */
  virtual ~HierarchicalAitkenPostProcessingTest () {};

  /**
   * @brief Empty.
   */
  virtual void setUp () {}

  /**
   * @brief Runs all tests.
   */
  virtual void run ();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;
};

}}} // namespace precice, cplscheme, tests

#endif /* PRECICE_CPLSCHEME_TESTS_HIERARCHICALAITKENPOSTPROCESSINGTEST_HPP_ */
