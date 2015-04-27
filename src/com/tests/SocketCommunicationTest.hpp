// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_SOCKETS

#ifndef PRECICE_COM_TESTS_SOCKETCOMMUNICATIONTEST_HPP_
#define PRECICE_COM_TESTS_SOCKETCOMMUNICATIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

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
  static tarch::logging::Log _log;

  void testSendAndReceive();

  void testParallelClient();

  void testReceiveFromAnyClient();

  //void testSendAndReceiveFromAnySender();
};

}}} // namespace precice, com, tests

#endif /* PRECICE_COM_TESTS_SOCKETCOMMUNICATIONTEST_HPP_ */

#endif // not PRECICE_NO_SOCKETS
