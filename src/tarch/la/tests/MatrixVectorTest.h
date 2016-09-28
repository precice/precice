#ifndef _LA_TARCH_TESTS_MATRIXVECTORTEST_H_
#define _LA_TARCH_TESTS_MATRIXVECTORTEST_H_

#include "tarch/tests/TestCase.h"

namespace tarch {
  namespace la {
    class MatrixVectorTest;
  }
}

/**
 * Provides tests for types Vector, DynamicVector and all Vector functionality.
 */
class tarch::la::MatrixVectorTest : public tarch::tests::TestCase
{
private:

  /**
   * Tests methods from MatrixVectorOperations.h.
   */
  void testMultiplication ();

  void testForwardSubstitution ();

  void testBackSubstitution ();

  void testSolveSystem3x3 ();

public:

  /**
   * Cosntructor.
   */
  MatrixVectorTest ();

  /**
   * Destructor, empty.
   */
  virtual ~MatrixVectorTest () {};

  /**
   * This routine is triggered by the TestCaseCollection
   */
  virtual void run();

  /**
   * Setup your test case.
   */
  virtual void setUp() {};
};

#endif /* _LA_TARCH_TESTS_MATRIXVECTORTEST_H_ */
