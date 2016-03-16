// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_TESTS_COMMUNICATIONMPIDIRECTTEST_HPP_
#define PRECICE_COM_TESTS_COMMUNICATIONMPIDIRECTTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace com {
namespace tests {

class MPIDirectCommunicationTest : public tarch::tests::TestCase
{
public:

  MPIDirectCommunicationTest();

  virtual ~MPIDirectCommunicationTest() {}

  /**
   * @brief Empty.
   */
  virtual void setUp() {}

  virtual void run();

private:

  static logging::Logger _log;

# ifndef PRECICE_NO_MPI

  void testSendReceiveTwoProcesses();

  void testSendReceiveThreeProcesses();

# endif // not PRECICE_NO_MPI
};

}}} // namespace precice, com, tests

#endif /* PRECICE_COM_TESTS_COMMUNICATIONMPIDIRECTTEST_HPP_ */
