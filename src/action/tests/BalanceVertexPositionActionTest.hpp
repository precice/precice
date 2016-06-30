#ifndef PRECICE_TESTS_BALANCEVERTEXPOSITIONACTIONTEST_HPP_
#define PRECICE_TESTS_BALANCEVERTEXPOSITIONACTIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace action {
namespace tests {

/**
 * @brief Provides tests for class precice::action::BalanceVertexPositionAction.
 */
class BalanceVertexPositionActionTest : public tarch::tests::TestCase
{
public:

  BalanceVertexPositionActionTest ();

  virtual ~BalanceVertexPositionActionTest () {}

  virtual void setUp() {}

  virtual void run ();

private:

  static logging::Logger _log;

  void testSmoothCircle ();

  void testSmoothSphere ();

  void testSmoothHexahedron ();

  void testConfiguration ();
};

}}} // namespace precice, action, tests

#endif /* PRECICE_TESTS_BALANCEVERTEXPOSITIONACTIONTEST_HPP_ */
