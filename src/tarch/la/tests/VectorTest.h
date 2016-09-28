#ifndef _LA_TARCH_TESTS_VECTORTEST_H_
#define _LA_TARCH_TESTS_VECTORTEST_H_

#include "tarch/tests/TestCase.h"

namespace tarch {
  namespace la {
    class VectorTest;
  }
}

/**
 * Provides tests for types Vector, DynamicVector and all Vector functionality.
 */
class tarch::la::VectorTest : public tarch::tests::TestCase
{
private:

  /**
   * Tests constructors.
   */
  void testConstruction();

  /**
   * Tests methods from VectorAssign.h and operator= from Vector types.
   */
  void testAssignment();

  /**
   * Tests methods from VectorOperations.h.
   */
  void testVectorOperations();

  /**
   * Tests methods from VectorScalarOperations.h.
   */
  void testVectorScalarOperations();

  /**
   * Tests methods from VectorVectorOperations.h.
   */
  void testVectorVectorOperations();

  /**
   * Tests wrapping a raw array with (static) vector semantics.
   */
  void testWrappedVector();

  /**
   * Tests methods from VectorCompare.h.
   */
  void testVectorCompare();

  void testVectorVectorCompare();

  void testVectorConversion();

  void testStdVector();

public:

  /**
   * Cosntructor.
   */
  VectorTest ();

  /**
   * Destructor, empty.
   */
  virtual ~VectorTest () {};

  /**
   * This routine is triggered by the TestCaseCollection
   */
  virtual void run();

  /**
   * Setup your test case.
   */
  virtual void setUp() {};
};

#endif /* _LA_TARCH_TESTS_VECTORTEST_H_ */
