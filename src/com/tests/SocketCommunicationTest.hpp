
#ifndef PRECICE_NO_SOCKETS

#ifndef PRECICE_COM_TESTS_SOCKETCOMMUNICATIONTEST_HPP_
#define PRECICE_COM_TESTS_SOCKETCOMMUNICATIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace com {
namespace tests {

/**
 * @brief Provides tests for class SocketCommunication.
 */
class SocketCommunicationTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  SocketCommunicationTest();

  /**
   * @brief Destructor, empty.
   */
  virtual ~SocketCommunicationTest() {}

  /**
   * @brief Empty.
   */
  virtual void setUp() {}

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  // @brief Logging device.
  static logging::Logger _log;

  void testSendAndReceive();

  void testParallelClient();

  void testReceiveFromAnyClient();

  //void testSendAndReceiveFromAnySender();
};

}}} // namespace precice, com, tests

#endif /* PRECICE_COM_TESTS_SOCKETCOMMUNICATIONTEST_HPP_ */

#endif // not PRECICE_NO_SOCKETS
