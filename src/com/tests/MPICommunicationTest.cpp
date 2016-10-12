#include "MPICommunicationTest.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "math/math.hpp"
#include "tarch/tests/TestCaseFactory.h"

using precice::utils::Parallel;

registerTest(precice::com::tests::MPICommunicationTest);

namespace precice {
namespace com {
namespace tests {

logging::Logger MPICommunicationTest::_log(
    "precice::com::tests::MPICommunicationTest");

MPICommunicationTest::MPICommunicationTest()
    : TestCase("precice::com::tests::MPICommunicationTest") {
}

void
MPICommunicationTest::run() {
# ifndef PRECICE_NO_MPI

  preciceTrace("run");

  Parallel::synchronizeProcesses();

  if (Parallel::getCommunicatorSize() > 1) {
    auto communicator = Parallel::getRestrictedCommunicator({0, 1});

    if (Parallel::getProcessRank() < 2) {
      Parallel::setGlobalCommunicator(communicator);
      testMethod(testSendAndReceiveString);
      testMethod(testSendAndReceiveVector);
      testMethod(testSendAndReceiveInteger);
      Parallel::setGlobalCommunicator(Parallel::getCommunicatorWorld());
    }
  }
# endif
}

#ifndef PRECICE_NO_MPI

void
MPICommunicationTest::testSendAndReceiveString() {
  preciceTrace("testSendAndReceiveString()");
  utils::Parallel::synchronizeProcesses();
  MPIPortsCommunication com;
  if (utils::Parallel::getProcessRank() == 0) {
    com.acceptConnection("A", "R", 0, 1);
    std::string msg("testOne");
    com.send(msg, 0);
    com.receive(msg, 0);
    validate(msg == std::string("testTwo"));

  } else if (utils::Parallel::getProcessRank() == 1) {
    com.requestConnection("A", "R", 0, 1);
    std::string msg;
    com.receive(msg, 0);
    validate(msg == std::string("testOne"));
    msg = "testTwo";
    com.send(msg, 0);
  }
  utils::Parallel::clearGroups();
}

void
MPICommunicationTest::testSendAndReceiveVector() {
  preciceTrace("testSendAndReceiveVector()");
  utils::Parallel::synchronizeProcesses();

  MPIPortsCommunication com;
  if (utils::Parallel::getProcessRank() == 0) {
    com.acceptConnection("A", "R", 0, 1);
    Eigen::Vector3d msg(1, 1, 1);
    com.send(msg.data(), msg.size(), 0);
    com.receive(msg.data(), msg.size(), 0);
    validate(math::equals(msg, Eigen::Vector3d::Constant(2)));
  } else if (utils::Parallel::getProcessRank() == 1) {
    com.requestConnection("A", "R", 0, 1);
    Eigen::Vector3d msg(0, 0, 0);
    com.receive(msg.data(), msg.size(), 0);
    validate(math::equals(msg, Eigen::Vector3d::Constant(1)));
    msg = Eigen::Vector3d::Constant(2);
    com.send(msg.data(), msg.size(), 0);
  }
  utils::Parallel::clearGroups();
}

void
MPICommunicationTest::testSendAndReceiveInteger() {
  preciceTrace("testSendAndReceiveInteger()");
  utils::Parallel::synchronizeProcesses();

  MPIPortsCommunication com;
  if (utils::Parallel::getProcessRank() == 0) {
    com.acceptConnection("A", "R", 0, 1);
    int msg = 1;
    com.send(msg, 0);
    com.receive(msg, 0);
    validate(msg == 2);
  } else if (utils::Parallel::getProcessRank() == 1) {
    com.requestConnection("A", "R", 0, 1);
    int msg = 0;
    com.receive(msg, 0);
    validate(msg == 1);
    msg = 2;
    com.send(msg, 0);
  }
  utils::Parallel::clearGroups();
}


#endif // not PRECICE_NO_MPI
}
}
} // namespace precice, com, tests
