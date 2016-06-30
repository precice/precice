#ifndef _LA_TARCH_TESTS_QRFACTORIZATIONTEST_H_
#define _LA_TARCH_TESTS_QRFACTORIZATIONTEST_H_

#include "tarch/tests/TestCase.h"
#include <tarch/la/DynamicMatrix.h>
#include <Eigen/Dense>

namespace precice {
  namespace cplscheme {
    namespace tests {
    class QRFactorizationTest;
    }
  }
}

namespace precice {
namespace cplscheme {
namespace tests {
/**
 * Provides tests for types Vector, DynamicVector and all Vector functionality.
 */
class QRFactorizationTest : public tarch::tests::TestCase
{
private:

  /**
   * Tests constructors.
   */
  void testQRFactorization ();
  
  void testQTQequalsIdentity(Eigen::MatrixXd& Q);
  void testQTQequalsIdentity(tarch::la::DynamicMatrix<double>& dynQ);
  
  void testQRequalsA(Eigen::MatrixXd& Q, Eigen::MatrixXd& R, Eigen::MatrixXd& A);

public:

  /**
   * Cosntructor.
   */
  QRFactorizationTest ();

  /**
   * Destructor, empty.
   */
  virtual ~QRFactorizationTest () {}

  /**
   * This routine is triggered by the TestCaseCollection
   */
  virtual void run();

  /**
   * Setup your test case.
   */
  virtual void setUp() {};
};

}}} // namespace precice, cplscheme, tests

#endif /* _LA_TARCH_TESTS_QRFACTORIZATIONTEST_H_ */
