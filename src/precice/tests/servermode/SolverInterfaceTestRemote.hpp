#ifndef PRECICE_ITESTS_COUPLINGINTERFACETESTREMOTE_HPP_
#define PRECICE_ITESTS_COUPLINGINTERFACETESTREMOTE_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include <string>
#include "precice/impl/SharedPointer.hpp"
#include "precice/SolverInterface.hpp"

namespace precice {
namespace tests {

class SolverInterfaceTestRemote : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  SolverInterfaceTestRemote();

  /**
   * @brief Destructor.
   */
  virtual ~SolverInterfaceTestRemote();

  /**
   * @brief Retrieves path to test directory.
   */
  virtual void setUp();

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  // @brief Logging device.
  static logging::Logger _log;

  // @brief Path to this directory.
  std::string _pathToTests;

  /**
   * @brief As SolverInterface::configure(), but without changing logging config.
   */
  void configureSolverInterface (
    const std::string& configFilename,
    SolverInterface&   interface );

  /**
   * @brief Runs two solver interfaces in coupling mode, one using a server.
   */
  void testCouplingModeWithOneServer();

  /**
   * @brief Two solvers in coupling mode, one in parallel using a server.
   */
  void testCouplingModeParallelWithOneServer();
};

}} // namespace precice, tests

#endif /* PRECICE_ITESTS_COUPLINGINTERFACETESTREMOTE_HPP_ */
