// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_TESTS_SCALARTEST_H_
#define _TARCH_LA_TESTS_SCALARTEST_H_

#include "tarch/tests/TestCase.h"

namespace tarch {
  namespace la {
    class ScalarTest;
  }
}

class tarch::la::ScalarTest : public tarch::tests::TestCase
{
private:

  void testComparison ();

  void testAbs ();

public:

  /**
   * Constructor.
   */
  ScalarTest ();

  /**
   * Destructor, empty.
   */
  virtual ~ScalarTest () {};

  /**
   * Runs all tests.
   */
  virtual void run ();

  /**
   * Sets up test environment.
   */
  virtual void setUp () {}
};

#endif /* _TARCH_LA_TESTS_SCALARTEST_H_ */
