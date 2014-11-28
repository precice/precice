// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NO_MPI

#include "MPICommunication.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace com {

tarch::logging::Log MPICommunication::
    _log ("precice::com::MPICommunication");


MPICommunication:: MPICommunication
(
  const MPI_Comm& communicator )
:
  Communication (),
  _communicator ( communicator )
{
  _rankOffset = 0;
}

void MPICommunication:: send
(
  const std::string& itemToSend,
  int                rankReceiver )
{
  preciceTrace2 ( "send()", itemToSend, rankReceiver );
  rankReceiver = rankReceiver - _rankOffset;
  assertion ( rankReceiver != ANY_SENDER );
  int length = itemToSend.size();
  preciceDebug ( "Message length: " << length );
  MPI_Send ( &length, 1, MPI_INT, rankReceiver, 0, _communicator );
  char * cstringMessage = new char[length+1];
  cstringMessage = const_cast<char*>(itemToSend.c_str());
  preciceDebug ( "Message: " + std::string(cstringMessage) );
  MPI_Send (cstringMessage, length+1, MPI_CHAR, rankReceiver, 0, _communicator);
}

void MPICommunication:: send
(
 int* itemsToSend,
 int  size,
 int  rankReceiver )
{
  preciceTrace1 ( "send(int*)", size );
  rankReceiver = rankReceiver - _rankOffset;
  assertion ( rankReceiver != ANY_SENDER );
  MPI_Send ( itemsToSend, size, MPI_INT, rankReceiver, 0, _communicator );
}

void MPICommunication:: send
(
 double* itemsToSend,
 int     size,
 int     rankReceiver )
{
  preciceTrace1 ( "send(double*)", size );
  rankReceiver = rankReceiver - _rankOffset;
  assertion ( rankReceiver != ANY_SENDER );
  MPI_Send ( itemsToSend, size, MPI_DOUBLE, rankReceiver, 0, _communicator );
}

void MPICommunication:: send
(
   double itemToSend,
   int    rankReceiver )
{
  preciceTrace2 ( "send(double)", itemToSend, rankReceiver );
  rankReceiver = rankReceiver - _rankOffset;
  assertion ( rankReceiver != ANY_SENDER );
  MPI_Send (&itemToSend, 1, MPI_DOUBLE, rankReceiver, 0, _communicator);
}

void MPICommunication:: send
(
   int itemToSend,
   int rankReceiver )
{
  preciceTrace2 ( "send(int)", itemToSend, rankReceiver );
  rankReceiver = rankReceiver - _rankOffset;
  assertion ( rankReceiver != ANY_SENDER );
  MPI_Send (&itemToSend, 1, MPI_INT, rankReceiver, 0, _communicator);
}

void MPICommunication:: send
(
  bool itemToSend,
  int  rankReceiver )
{
  preciceTrace2 ( "send(bool)", itemToSend, rankReceiver );
  rankReceiver = rankReceiver - _rankOffset;
  assertion ( rankReceiver != ANY_SENDER );
  int buffer = itemToSend ? 1 : 0;
  MPI_Send ( &buffer, 1, MPI_INT, rankReceiver, 0, _communicator );
}

int MPICommunication:: receive
(
  std::string& itemToReceive,
  int          rankSender )
{
  preciceTrace2 ( "receive(string)", itemToReceive, rankSender );
  rankSender = rankSender - _rankOffset;
  int length;
  MPI_Status status1;
  rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
  MPI_Recv (&length, 1, MPI_INT, rankSender, 0, _communicator, &status1);
  rankSender = status1.MPI_SOURCE;
  preciceDebug ( "Stringlength = " << length );
  char * cstringMessage = new char[length+1];
  MPI_Status status2;
  MPI_Recv (cstringMessage, length+1, MPI_CHAR, rankSender, 0,
            _communicator, &status2);
  itemToReceive = std::string (cstringMessage);
  preciceDebug ( "Received \"" << itemToReceive << "\" from rank " << rankSender );
  return rankSender;
}

int MPICommunication:: receive
(
  int* itemsToReceive,
  int  size,
  int  rankSender )
{
  preciceTrace1 ( "receive(int*)", size );
  rankSender = rankSender - _rankOffset;
  rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
  MPI_Status status;
  MPI_Recv ( itemsToReceive, size, MPI_INT, rankSender, 0, _communicator, &status );
  return status.MPI_SOURCE;
}

int MPICommunication:: receive
(
  double * itemsToReceive,
  int      size,
  int      rankSender )
{
  preciceTrace1 ( "receive(double*)", size );
  rankSender = rankSender - _rankOffset;
  rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
  MPI_Status status;
  MPI_Recv ( itemsToReceive, size, MPI_DOUBLE, rankSender, 0, _communicator, &status );
  return status.MPI_SOURCE;
}

int MPICommunication:: receive
(
   double & itemToReceive,
   int      rankSender )
{
   preciceTrace1 ( "receive(double)", rankSender );
   rankSender = rankSender - _rankOffset;
   rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
   MPI_Status status;
   MPI_Recv (&itemToReceive, 1, MPI_DOUBLE, rankSender, 0, _communicator, &status);
   preciceDebug ( "Received " << itemToReceive << " from rank " << status.MPI_SOURCE );
   return status.MPI_SOURCE;
}

int MPICommunication:: receive
(
   int & itemToReceive,
   int   rankSender )
{
   preciceTrace1 ( "receive(int)", rankSender );
   rankSender = rankSender - _rankOffset;
   rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
   MPI_Status status;
   MPI_Recv (&itemToReceive, 1, MPI_INT, rankSender, 0,
             _communicator, &status);
   preciceDebug ( "Received " << itemToReceive << " from rank " << status.MPI_SOURCE );
   return status.MPI_SOURCE;
}

int MPICommunication:: receive
(
   bool & itemToReceive,
   int    rankSender )
{
   preciceTrace1 ( "receive(bool)", rankSender );
   rankSender = rankSender - _rankOffset;
   rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
   MPI_Status status;
   int buffer = -1;
   MPI_Recv ( &buffer, 1, MPI_INT, rankSender, 0, _communicator, &status );
   assertion ( buffer != -1 );
   itemToReceive = (buffer == 1) ? true : false;
   preciceDebug ( "Received " << itemToReceive << " from rank " << status.MPI_SOURCE );
   return status.MPI_SOURCE;
}

}} // namespace precice, com

#endif // not PRECICE_NO_MPI
