// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_ACTION_TESTS_PYTHONACTIONTEST_HPP_
#define PRECICE_ACTION_TESTS_PYTHONACTIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace action {
namespace tests {

/**
 * @brief Provides tests for class PythonAction.
 */
class PythonActionTest : public tarch::tests::TestCase
{
public:

  PythonActionTest ();

  virtual ~PythonActionTest () {}

  virtual void setUp () {}

  virtual void run ();

private:

  static tarch::logging::Log _log;

  void testAllMethods();

  void testOmitMethods();
};

}}} // namespace precice, action, tests

#endif /* PRECICE_ACTION_TESTS_PYTHONACTIONTEST_HPP_ */
