// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_SOCKETS

#include "SocketCommunicationTest.hpp"
#include "com/SocketCommunication.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/DynamicVector.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::com::tests::SocketCommunicationTest)

namespace precice {
namespace com {
namespace tests {

tarch::logging::Log SocketCommunicationTest:: _log ("precice::com::tests::SocketCommunicationTest");

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
    if ( utils::Parallel::getProcessRank() < 3 ){
      testMethod ( testReceiveFromAnyClient );
    }
    utils::Parallel::synchronizeProcesses(); // Necessary for sockets
  }
}

void SocketCommunicationTest:: testSendAndReceive()
{
  preciceTrace ( "testSendAndReceiveString()" );
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
      utils::DynVector msg(3, 0.0);
      com.receive (tarch::la::raw(msg), msg.size(), 0);
      validate ( tarch::la::equals(msg, utils::Vector3D(1.0)) );
      msg = utils::Vector3D (2.0);
      com.send (tarch::la::raw(msg), msg.size(), 0);
    }
    {
      tarch::la::DynamicVector<int> msg(4, 0);
      com.receive (tarch::la::raw(msg), msg.size(), 0);
      validate ( tarch::la::equals(msg, tarch::la::DynamicVector<int>(4, 1)) );
      assign(msg) = 0;
      com.send (tarch::la::raw(msg), msg.size(), 0);
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
      utils::DynVector msg(3, 1.0);
      com.send (tarch::la::raw(msg), msg.size(), 0);
      com.receive (tarch::la::raw(msg), msg.size(), 0);
      validate ( tarch::la::equals(msg, utils::Vector3D(2.0)) );
    }
    {
      tarch::la::DynamicVector<int> msg(4, 1);
      com.send (tarch::la::raw(msg), msg.size(), 0);
      com.receive (tarch::la::raw(msg), msg.size(), 0);
      validate ( tarch::la::equals(msg, tarch::la::DynamicVector<int>(4, 0)) );
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
  preciceTrace ( "testParallelClient()" );
  SocketCommunication com;
  int rank = utils::Parallel::getProcessRank();
  if ( rank == 0 ){
    preciceDebug("branch rank 0");
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
    preciceDebug("branch rank 1, 2");
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
      assertion1 ( rank == 2, rank );
      com.receive ( receiveMsg, 0 );
      validateEquals ( receiveMsg, 2 );
    }
    com.closeConnection();
  }
}

void SocketCommunicationTest:: testReceiveFromAnyClient()
{
  preciceTrace ( "testReceiveFromAnyClient()" );
  SocketCommunication com;
  int rank = utils::Parallel::getProcessRank();
  int rank0 = 0;
  int rank1 = 1;
  int rank2 = 0;
  if ( rank == 0 ){
    com.requestConnection("server", "client", 0, 2);
    int msg = 0;
    com.send(msg, rank2);
    int sender = com.receive(msg, rank2);
    validateEquals(sender, rank2);
    validateEquals(msg, 1);
    msg = 10;
    com.send(msg, rank2);
  }
  else if ( rank == 1 ){
    com.requestConnection("server", "client", 1, 2);
    int msg = -1;
    int sender = com.receive(msg, rank2);
    validateEquals(sender, rank2);
    validateEquals(msg, 2);
    msg = 20;
    com.send(msg, rank0);
  }
  else if ( rank == 2 ){
    com.acceptConnection("server", "client", 0, 1);
    int msg = -1;
    int sender = com.receive(msg, Communication::ANY_SENDER);
    validateEquals(sender, rank0);
    validateEquals(msg, 0);
    msg = 2;
    com.send(msg, rank1);
    sender = com.receive(msg, Communication::ANY_SENDER);
    validateEquals(sender, rank1);
    validateEquals(msg, 20);
    msg = 1;
    com.send(msg, rank0);
    sender = com.receive(msg, Communication::ANY_SENDER);
    validateEquals(sender, rank0);
    validateEquals(msg, 10);
  }
}

}}} // namespace precice, com, tests

#endif // not PRECICE_NO_SOCKETS
