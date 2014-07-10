// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _LA_TARCH_TESTS_GRAMSCHMIDTTEST_H_
#define _LA_TARCH_TESTS_GRAMSCHMIDTTEST_H_

#include "tarch/tests/TestCase.h"

namespace tarch {
  namespace la {
    class GramSchmidtTest;
  }
}

/**
 * Provides tests for types Vector, DynamicVector and all Vector functionality.
 */
class tarch::la::GramSchmidtTest : public tarch::tests::TestCase
{
private:

  /**
   * Tests constructors.
   */
  void testModifiedGramSchmidt ();

public:

  /**
   * Cosntructor.
   */
  GramSchmidtTest ();

  /**
   * Destructor, empty.
   */
  virtual ~GramSchmidtTest () {}

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
