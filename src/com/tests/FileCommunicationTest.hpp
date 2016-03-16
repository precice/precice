// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_TESTS_COMMUNICATIONFILETEST_HPP_
#define PRECICE_COM_TESTS_COMMUNICATIONFILETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

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
  static logging::Logger _log;

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
