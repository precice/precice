#pragma once

#include "tarch/tests/TestCase.h"
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

class QRFactorizationTest : public tarch::tests::TestCase
{
private:

  /**
   * Tests constructors.
   */
  void testQRFactorization ();
  
  void testQTQequalsIdentity(Eigen::MatrixXd& Q);

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

