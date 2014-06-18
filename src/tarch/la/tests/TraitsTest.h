// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_TESTS_TRAITSTEST_H_
#define _TARCH_LA_TESTS_TRAITSTEST_H_

#include "tarch/tests/TestCase.h"

namespace tarch {
  namespace la {
    class TraitsTest;
  }
}

class tarch::la::TraitsTest : public tarch::tests::TestCase
{
private:

//  void testIsEqual ();

  void testEqualScalars ();

public:

  TraitsTest ();

  virtual ~TraitsTest () {};

  /**
   * This routine is triggered by the TestCaseCollection
   */
  virtual void run();

  /**
   * Setup your test case.
   */
  virtual void setUp() {};
};

#endif /* _TARCH_LA_TESTS_TRAITSTEST_H_ */
