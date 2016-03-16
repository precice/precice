// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_TESTS_ABSOLUTECONVERGENCEMEASURETEST_HPP_
#define PRECICE_CPLSCHEME_TESTS_ABSOLUTECONVERGENCEMEASURETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace cplscheme {
namespace tests {

class AbsoluteConvergenceMeasureTest : public tarch::tests::TestCase
{
public:

  AbsoluteConvergenceMeasureTest ();

  virtual ~AbsoluteConvergenceMeasureTest () {};

  /**
   * @brief Empty.
   */
  virtual void setUp () {}

  virtual void run ();

private:

  static logging::Logger _log;

  //   void testMeasureVectorData ();

  void testMeasureData ();

  //   void testMeasureIntegerData ();
};

}}} // namespace precice, cplscheme, tests

#endif // PRECICE_CPLSCHEME_TESTS_ABSOLUTECONVERGENCEMEASURETEST_HPP_
