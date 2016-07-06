#ifndef PRECICE_CPLSCHEME_TESTS_UNCOUPLEDCOUPLINGSCHEMETEST_HPP_
#define PRECICE_CPLSCHEME_TESTS_UNCOUPLEDCOUPLINGSCHEMETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace cplscheme {
namespace tests {

/**
 * @brief Provides tests for class UncoupledScheme.
 */
class UncoupledSchemeTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  UncoupledSchemeTest ();

  /**
   * @brief Destructor.
   */
  virtual ~UncoupledSchemeTest () {}

  /**
   * @brief Prepares running test cases, empty.
   */
  virtual void setUp () {}

  /**
   * @brief Runs all tests.
   */
  virtual void run ();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  /**
   * @brief Only test case so far...
   */
  void testBasics ();
};

}}} // namespace precice, cplscheme, tests

#endif /* PRECICE_CPLSCHEME_TESTS_UNCOUPLEDCOUPLINGSCHEMETEST_HPP_ */
