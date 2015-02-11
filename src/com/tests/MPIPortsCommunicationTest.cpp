// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#include "MPIPortsCommunicationTest.hpp"

#include "com/MPIPortsCommunication.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"
#include "tarch/tests/TestCaseFactory.h"

using precice::utils::Parallel;

registerTest(precice::com::tests::MPIPortsCommunicationTest);

namespace precice {
namespace com {
namespace tests {

tarch::logging::Log MPIPortsCommunicationTest::_log(
    "precice::com::tests::MPIPortsCommunicationTest");

MPIPortsCommunicationTest::MPIPortsCommunicationTest()
    : TestCase("precice::com::tests::MPIPortsCommunicationTest") {
}

void
MPIPortsCommunicationTest::run() {
  preciceTrace("run");

  Parallel::synchronizeProcesses();

  if (Parallel::getCommunicatorSize() > 1) {
    auto communicator = Parallel::getRestrictedCommunicator({0, 1});

    if (Parallel::getProcessRank() < 2) {
      Parallel::setGlobalCommunicator(communicator);
      testMethod(testSendReceiveTwoProcesses);
      Parallel::setGlobalCommunicator(Parallel::getCommunicatorWorld());
    }
  }
}

void
MPIPortsCommunicationTest::testSendReceiveTwoProcesses() {
  preciceTrace("testSendReceiveTwoProcesses()");

  assertion(Parallel::getCommunicatorSize() == 2);

  validateEquals(Parallel::getCommunicatorSize(), 2);

  MPIPortsCommunication communication;

  std::string nameEven("even");
  std::string nameOdd("odd");

  switch (Parallel::getProcessRank()) {
  case 0: {
    communication.acceptConnection(nameEven, nameOdd, 0, 1);

    int message = 1;

    communication.send(message, 0);
    communication.receive(message, 0);

    validateEquals(message, 2);

    communication.closeConnection();

    break;
  }
  case 1: {
    communication.requestConnection(nameEven, nameOdd, 0, 1);

    int message = -1;

    communication.receive(message, 0);

    validateEquals(message, 1);

    message = 2;

    communication.send(message, 0);

    communication.closeConnection();

    break;
  }
  }

  Parallel::synchronizeProcesses();
}
}
}
} // namespace precice, com, tests

#endif // not PRECICE_NO_MPI
