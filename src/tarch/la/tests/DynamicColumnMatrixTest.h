#ifndef _LA_TARCH_TESTS_DYNAMICCOLUMNMATRIXTEST_H_
#define _LA_TARCH_TESTS_DYNAMICCOLUMNMATRIXTEST_H_

#include "tarch/tests/TestCase.h"

namespace tarch {
  namespace la {
    class DynamicColumnMatrixTest;
  }
}

/**
 * Provides tests for types Vector, DynamicVector and all Vector functionality.
 */
class tarch::la::DynamicColumnMatrixTest : public tarch::tests::TestCase
{
private:

  /**
   * Tests of basic matrix constructors and accessors.
   */
  void testBasics ();

  /**
   * Tests of column manipulation operations.
   */
  void testColumnManipulations ();

public:

  /**
   * Cosntructor.
   */
  DynamicColumnMatrixTest ();

  /**
   * Destructor, empty.
   */
  virtual ~DynamicColumnMatrixTest () {}

  /**
   * Empty.
   */
  virtual void setUp() {}

  /**
   * This routine is triggered by the TestCaseCollection
   */
  virtual void run();
};

#endif /* _LA_TARCH_TESTS_DYNAMICCOLUMNMATRIXTEST_H_ */
