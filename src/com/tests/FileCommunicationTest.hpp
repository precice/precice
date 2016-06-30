#ifndef PRECICE_COM_TESTS_COMMUNICATIONFILETEST_HPP_
#define PRECICE_COM_TESTS_COMMUNICATIONFILETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace com {
namespace tests {

/**
 * @brief Provides tests for class FileCommunication.
 */
class FileCommunicationTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  FileCommunicationTest ();

  /**
   * @brief Destructor, empty.
   */
  virtual ~FileCommunicationTest() {}

  /**
   * @brief Prepares execution of tests, empty.
   */
  virtual void setUp() {}

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  /**
   * @brief Tests a single send/receive of all possible message types.
   */
  void testSimpleSendReceive();

  /**
   * @brief Tests repeated send/receives of one message type.
   */
  void testMultipleExchanges();
};

}}} // namespace precice, com, tests

#endif /* PRECICE_COM_TESTS_COMMUNICATIONFILETEST_HPP_ */
