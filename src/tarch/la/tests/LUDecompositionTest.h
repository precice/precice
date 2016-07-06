#ifndef LUDECOMPOSITIONTEST_H_
#define LUDECOMPOSITIONTEST_H_

#include "tarch/tests/TestCase.h"

namespace tarch {
  namespace la {
    class LUDecompositionTest;
  }
}

class tarch::la::LUDecompositionTest : public tarch::tests::TestCase
{
private:

  void testLUNoPivoting();

  void testLU();

public:

  /**
   * Constructor.
   */
  LUDecompositionTest();

  /**
   * Destructor, empty.
   */
  virtual ~LUDecompositionTest() {}

  /**
   * This routine is triggered by the TestCaseCollection
   */
  virtual void run();

  /**
   * Setup your test case.
   */
  virtual void setUp() {}
};

#endif /* LUDECOMPOSITIONTEST_H_ */
