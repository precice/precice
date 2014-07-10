// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _LA_TARCH_TESTS_MATRIXTEST_H_
#define _LA_TARCH_TESTS_MATRIXTEST_H_

#include "tarch/tests/TestCase.h"

namespace tarch {
  namespace la {
    class MatrixTest;
  }
}

/**
 * Provides tests for types Vector, DynamicVector and all Vector functionality.
 */
class tarch::la::MatrixTest : public tarch::tests::TestCase
{
private:

  /**
   * Tests constructors.
   */
  void testConstruction ();

  /**
   * Tests methods from MatrixAssign.h.
   */
  void testAssignment ();

  /**
   * Tests methods from MatrixOperations.h.
   */
  void testMatrixOperations ();

  /**
   * Tests methods from MatrixMatrixOperations.h.
   */
  void testMatrixMatrixOperations ();

  /**
   * Tests operations with TransposedMatrix.h.
   */
  void testTransposedMatrix ();

public:

  /**
   * Cosntructor.
   */
  MatrixTest ();

  /**
   * Destructor, empty.
   */
  virtual ~MatrixTest () {};

  /**
   * This routine is triggered by the TestCaseCollection
   */
  virtual void run();

  /**
   * Setup your test case.
   */
  virtual void setUp() {};
};

#endif /* _LA_TARCH_TESTS_MATRIXTEST_H_ */
