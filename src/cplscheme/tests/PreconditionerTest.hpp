#pragma once

#include "tarch/tests/TestCase.h"
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
  typedef Eigen::VectorXd DataValues;
    
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

  void testMultilpleMeshes ();

  DataValues _data;
  DataValues _res;
  DataValues _compareDataRes;
  DataValues _compareDataResSum;
  DataValues _compareDataResSum2;
  DataValues _compareDataValue;
  DataValues _compareDataConstant;

  static logging::Logger _log;
};

}}} // namespace precice, cplscheme, tests
