// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_TESTS_BALANCEVERTEXPOSITIONACTIONTEST_HPP_
#define PRECICE_TESTS_BALANCEVERTEXPOSITIONACTIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace action {
namespace tests {

/**
 * @brief Provides tests for class precice::action::BalanceVertexPositionAction.
 */
class BalanceVertexPositionActionTest : public tarch::tests::TestCase
{
public:

  BalanceVertexPositionActionTest ();

  virtual ~BalanceVertexPositionActionTest () {}

  virtual void setUp() {}

  virtual void run ();

private:

  static logging::Logger _log;

  void testSmoothCircle ();

  void testSmoothSphere ();

  void testSmoothHexahedron ();

  void testConfiguration ();
};

}}} // namespace precice, action, tests

#endif /* PRECICE_TESTS_BALANCEVERTEXPOSITIONACTIONTEST_HPP_ */
