// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "SimpleCommunication.hpp"
#include "com/Communication.hpp"

namespace precice {
namespace m2n {

tarch::logging::Log SimpleCommunication:: _log("precice::m2n::SimpleCommunication");

SimpleCommunication:: SimpleCommunication
(
  com::PtrCommunication com )
:
  _com(com)
{}

SimpleCommunication:: ~SimpleCommunication()
{
  if (isConnected()){
    closeConnection();
  }
}

bool SimpleCommunication:: isConnected()
{
  return _com->isConnected();
}

void SimpleCommunication:: acceptConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                acceptorProcessRank,
  int                acceptorCommunicatorSize )
{
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);
  _com->acceptConnection(nameAcceptor, nameRequester, acceptorProcessRank, acceptorCommunicatorSize);
}

void SimpleCommunication:: requestConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                requesterProcessRank,
  int                requesterCommunicatorSize )
{
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);
  _com->requestConnection(nameAcceptor, nameRequester, requesterProcessRank, requesterCommunicatorSize);
}

void SimpleCommunication:: closeConnection()
{
  preciceTrace("closeConnection()");
  _com->closeConnection();
}

com::PtrCommunication SimpleCommunication:: getMasterCommunication()
{
  return _com;
}

void SimpleCommunication:: startSendPackage
(
  int rankReceiver )
{
  _com->startSendPackage(rankReceiver);
}

void SimpleCommunication:: finishSendPackage()
{
  _com->finishSendPackage();
}

int SimpleCommunication:: startReceivePackage
(
  int rankSender )
{
  preciceTrace1("startReceivePackage()", rankSender);
  return _com->startReceivePackage(rankSender);
}

void SimpleCommunication:: finishReceivePackage()
{
  _com->finishReceivePackage();
}

void SimpleCommunication:: sendMaster
(
  const std::string& itemToSend,
  int                rankReceiver )
{
  preciceTrace2("send(string)", itemToSend, rankReceiver);
  _com->send(itemToSend, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  int* itemsToSend,
  int  size,
  int  rankReceiver )
{
  preciceTrace2("send(int*)", size, rankReceiver);
  _com->send(itemsToSend, size, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  double* itemsToSend,
  int     size,
  int     rankReceiver )
{
  preciceTrace2("send(double*)", size, rankReceiver);
  _com->send(itemsToSend, size, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  double itemToSend,
  int    rankReceiver )
{
  preciceTrace2("send(double)", itemToSend, rankReceiver);
  _com->send(itemToSend, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  int itemToSend,
  int rankReceiver )
{
  preciceTrace2("send(int)", itemToSend, rankReceiver);
  _com->send(itemToSend, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  bool itemToSend,
  int  rankReceiver )
{
  preciceTrace2("send(bool)", itemToSend, rankReceiver);
  _com->send(itemToSend, rankReceiver);
}

int SimpleCommunication:: receiveMaster
(
  std::string& itemToReceive,
  int          rankSender )
{
  preciceTrace1("receive(string)", rankSender);
  return _com->receive(itemToReceive, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  int* itemsToReceive,
  int  size,
  int  rankSender )
{
  preciceTrace2("receive(int*)", size, rankSender);
  return _com->receive(itemsToReceive, size, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  double* itemsToReceive,
  int     size,
  int     rankSender )
{
  preciceTrace2("receive(double*)", size, rankSender);
  return _com->receive(itemsToReceive, size, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  double& itemToReceive,
  int     rankSender )
{
  preciceTrace1("receive(double)", rankSender);
  return _com->receive(itemToReceive, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  int& itemToReceive,
  int  rankSender )
{
  preciceTrace1("receive(int)", rankSender);
  return _com->receive(itemToReceive, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  bool& itemToReceive,
  int   rankSender )
{
  preciceTrace1("receive(bool)", rankSender);
  return _com->receive(itemToReceive, rankSender);
}


void SimpleCommunication:: sendAll (
  double* itemsToSend,
  int     size,
  int     rankReceiver )
{
  _com->send(itemsToSend, size, rankReceiver);
}


/**
 * @brief Receives an array of double values.
 *
 * @return Rank of sender, which is useful when ANY_SENDER is used.
 */
int SimpleCommunication:: receiveAll (
  double* itemsToReceive,
  int     size,
  int     rankSender )
{
  return _com->receive(itemsToReceive, size, rankSender);
}

}} // namespace precice, m2n

