// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_DIMENSIONSTEST_HPP_
#define PRECICE_UTILS_DIMENSIONSTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace utils {
namespace tests {

/**
 * Provides tests for methods/classes in file utils/Dimensions.
 */
class DimensionsTest : public tarch::tests::TestCase
{
public:

  DimensionsTest();

  virtual ~DimensionsTest() {}

  virtual void setUp() {}

  virtual void run();

private:

  static tarch::logging::Log _log;

  void testLinearizeDelinearize();

//  void testGetHyperfaceCornerIndices();
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_DIMENSIONSTEST_HPP_ */
