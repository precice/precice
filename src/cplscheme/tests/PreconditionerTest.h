// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#pragma once

#include "tarch/tests/TestCase.h"
#include <tarch/la/DynamicColumnMatrix.h>
#include <tarch/la/DynamicVector.h>
#include <Eigen/Dense>
#include "utils/Globals.hpp"

namespace precice {
  namespace cplscheme {
    namespace tests {
    class PreconditionerTest;
    }
  }
}

namespace precice {
namespace cplscheme {
namespace tests {
/**
 * Provides tests for all preconditioner variants
 */
class PreconditionerTest : public tarch::tests::TestCase
{
public:

  //see post-processing definitions
  typedef tarch::la::DynamicVector<double> DataValues;
  typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;
  typedef Eigen::MatrixXd EigenMatrix;

  /**
   * Constructor.
   */
  PreconditionerTest ();

  /**
   * Destructor, empty.
   */
  virtual ~PreconditionerTest () {}

  /**
   * This routine is triggered by the TestCaseCollection
   */
  virtual void run();

  /**
   * Setup your test case.
   */
  virtual void setUp();

private:

  /**
   * Tests residual preconditioner.
   */
  void testResPreconditioner ();

  /**
   * Tests residual sum preconditioner.
   */
  void testResSumPreconditioner ();

  /**
   * Tests value preconditioner.
   */
  void testValuePreconditioner ();

  /**
   * Tests constant preconditioner.
   */
  void testConstPreconditioner ();


  /**
   * Test preconditioner on squared matrix (aka old Jacobian) as well as tall and skinny matrix (aka V or W)
   */
  void testParallelMatrixScaling ();

  void validateVector(DataValues& data, DataValues& compare);

  DataValues _data;
  DataValues _res;
  DataValues _compareDataRes;
  DataValues _compareDataResSum;
  DataValues _compareDataValue;
  DataValues _compareDataConstant;

  static tarch::logging::Log _log;
};

}}} // namespace precice, cplscheme, tests
