
#ifndef PRECICE_NO_SOCKETS

#include "SocketCommunicationTest.hpp"
#include "com/SocketCommunication.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "math/math.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::com::tests::SocketCommunicationTest)

namespace precice {
namespace com {
namespace tests {

logging::Logger SocketCommunicationTest:: _log("precice::com::tests::SocketCommunicationTest");

SocketCommunicationTest:: SocketCommunicationTest()
:
  TestCase ("precice::com::tests::SocketCommunicationTest")
{}

void SocketCommunicationTest:: run()
{
  if ( utils::Parallel::getCommunicatorSize() >= 2 ){
    if ( utils::Parallel::getProcessRank() < 2 ){
      testMethod ( testSendAndReceive );
    }
    utils::Parallel::synchronizeProcesses(); // Necessary for sockets
  }
  if ( utils::Parallel::getCommunicatorSize() >= 3 ){
    if ( utils::Parallel::getProcessRank() < 3 ){
      testMethod ( testParallelClient );
    }
    utils::Parallel::synchronizeProcesses(); // Necessary for sockets
  }
}

void SocketCommunicationTest:: testSendAndReceive()
{
  TRACE();
  SocketCommunication com;
  if ( utils::Parallel::getProcessRank() == 0 ){
    com.acceptConnection("process0", "process1", 0, 1);
    {
      std::string msg ("testOne");
      com.send (msg, 0);
      com.receive (msg, 0);
      validate ( msg == std::string("testTwo") );
    }
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(0);
      com.receive (msg.data(), msg.size(), 0);
      validate ( math::equals(msg, Eigen::Vector3d::Constant(1)) );
      msg = Eigen::Vector3d::Constant(2);
      com.send (msg.data(), msg.size(), 0);
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(0);
      com.receive (msg.data(), msg.size(), 0);
      validate ( math::equals(msg, Eigen::Vector4i::Constant(1)) );
      msg = Eigen::Vector4i::Constant(0);
      com.send (msg.data(), msg.size(), 0);
    }
    {
      double msg = 0.0;
      com.send (msg, 0);
      com.receive (msg, 0);
      validateNumericalEquals ( msg, 1.0 );
    }
    {
      int msg = 1;
      com.send (msg, 0);
      com.receive (msg, 0);
      validate ( msg == 2 );
    }
    {
      bool msg = true;
      com.send (msg, 0);
      com.receive (msg, 0);
      validate ( msg == false );
    }
    com.closeConnection();
  }
  else if ( utils::Parallel::getProcessRank() == 1 ){
    com.requestConnection("process0", "process1", 0, 1);
    {
      std::string msg;
      com.receive (msg, 0);
      validate ( msg == std::string("testOne") );
      msg = "testTwo";
      com.send (msg, 0);
    }
    {
      Eigen::Vector3d msg = Eigen::Vector3d::Constant(1);
      com.send (msg.data(), msg.size(), 0);
      com.receive (msg.data(), msg.size(), 0);
      validate ( math::equals(msg, Eigen::Vector3d::Constant(2)) );
    }
    {
      Eigen::Vector4i msg = Eigen::Vector4i::Constant(1);
      com.send (msg.data(), msg.size(), 0);
      com.receive (msg.data(), msg.size(), 0);
      validate ( math::equals(msg, Eigen::Vector4i::Zero()) );
    }
    {
      double msg = 1.0;
      com.receive (msg, 0);
      validateNumericalEquals ( msg, 0.0 );
      msg = 1.0;
      com.send (msg, 0);
    }
    {
      int msg = 0;
      com.receive (msg, 0);
      validate ( msg == 1 );
      msg = 2;
      com.send (msg, 0);
    }
    {
      bool msg = false;
      com.receive (msg, 0);
      validate ( msg == true );
      msg = false;
      com.send (msg, 0);
    }
    com.closeConnection();
  }
}

void SocketCommunicationTest:: testParallelClient()
{
  TRACE();
  SocketCommunication com;
  int rank = utils::Parallel::getProcessRank();
  if ( rank == 0 ){
    DEBUG("branch rank 0");
    com.acceptConnection("server", "client", 0, 1);
    validateEquals ( com.getRemoteCommunicatorSize(), 2 );
    std::string msg;
    com.receive (msg, 0);
    validate ( msg == std::string("process 0") );
    com.receive (msg, 1);
    validate ( msg == std::string("process 1") );
    int sendMsg = 1;
    com.send ( sendMsg, 0 );
    sendMsg = 2;
    com.send ( sendMsg, 1 );
    com.closeConnection();
  }
  else if ( (rank == 1) || (rank == 2) ){
    DEBUG("branch rank 1, 2");
    com.requestConnection("server", "client", rank-1, 2);
    validateEquals ( com.getRemoteCommunicatorSize(), 1 );
    std::ostringstream rankMsg;
    rankMsg << "process " << rank-1;
    com.send (rankMsg.str(), 0);
    int receiveMsg = 0;
    if ( rank == 1 ){
      com.receive ( receiveMsg, 0 );
      validateEquals ( receiveMsg, 1 );
    }
    else {
      assertion ( rank == 2, rank );
      com.receive ( receiveMsg, 0 );
      validateEquals ( receiveMsg, 2 );
    }
    com.closeConnection();
  }
}

}}} // namespace precice, com, tests

#endif // not PRECICE_NO_SOCKETS
