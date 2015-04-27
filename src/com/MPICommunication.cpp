// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#include "MPICommunication.hpp"

#include "MPIRequest.hpp"

#include "utils/Globals.hpp"

template <size_t>
struct MPI_Select_unsigned_integer_datatype;

template <>
struct MPI_Select_unsigned_integer_datatype<1> {
  static constexpr MPI_Datatype const datatype = MPI_UNSIGNED_CHAR;
};

template <>
struct MPI_Select_unsigned_integer_datatype<2> {
  static constexpr MPI_Datatype const datatype = MPI_UNSIGNED_SHORT;
};

template <>
struct MPI_Select_unsigned_integer_datatype<4> {
  static constexpr MPI_Datatype const datatype = MPI_UNSIGNED;
};

template <>
struct MPI_Select_unsigned_integer_datatype<8> {
  static constexpr MPI_Datatype const datatype = MPI_UNSIGNED_LONG;
};

#define MPI_BOOL MPI_Select_unsigned_integer_datatype<sizeof(bool)>::datatype

namespace precice {
namespace com {
tarch::logging::Log MPICommunication::_log("precice::com::MPICommunication");

MPICommunication::MPICommunication() {
}

void
MPICommunication::send(std::string const& itemToSend, int rankReceiver) {
  preciceTrace2("send()", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion(rankReceiver != ANY_SENDER);
  int length = itemToSend.size();
  preciceDebug("Message length: " << length);
  // const_cast is needed because MPI_Send expects a void* as first argument.
  char* cstr = const_cast<char*>(itemToSend.c_str());
  preciceDebug("Message: " + std::string(cstr));
  MPI_Send(cstr,
           length + 1,
           MPI_CHAR,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

void
MPICommunication::send(int* itemsToSend, int size, int rankReceiver) {
  preciceTrace1("send(int*)", size);
  rankReceiver = rankReceiver - _rankOffset;
  assertion(rankReceiver != ANY_SENDER);
  MPI_Send(itemsToSend,
           size,
           MPI_INT,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

Request::SharedPointer
MPICommunication::aSend(int* itemsToSend, int size, int rankReceiver) {
  preciceTrace1("aSend(int*)", size);
  rankReceiver = rankReceiver - _rankOffset;
  assertion(rankReceiver != ANY_SENDER);

  MPI_Request request;

  MPI_Isend(itemsToSend,
            size,
            MPI_INT,
            rank(rankReceiver),
            0,
            communicator(rankReceiver),
            &request);

  return Request::SharedPointer(new MPIRequest(request));
}

void
MPICommunication::send(double* itemsToSend, int size, int rankReceiver) {
  preciceTrace1("send(double*)", size);
  rankReceiver = rankReceiver - _rankOffset;
  assertion(rankReceiver != ANY_SENDER);
  MPI_Send(itemsToSend,
           size,
           MPI_DOUBLE,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

Request::SharedPointer
MPICommunication::aSend(double* itemsToSend, int size, int rankReceiver) {
  preciceTrace1("aSend(double*)", size);
  rankReceiver = rankReceiver - _rankOffset;
  assertion(rankReceiver != ANY_SENDER);

  MPI_Request request;

  MPI_Isend(itemsToSend,
            size,
            MPI_DOUBLE,
            rank(rankReceiver),
            0,
            communicator(rankReceiver),
            &request);

  return Request::SharedPointer(new MPIRequest(request));
}

void
MPICommunication::send(double itemToSend, int rankReceiver) {
  preciceTrace2("send(double)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion(rankReceiver != ANY_SENDER);
  MPI_Send(&itemToSend,
           1,
           MPI_DOUBLE,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

Request::SharedPointer
MPICommunication::aSend(double* itemToSend, int rankReceiver) {
  return aSend(itemToSend, 1, rankReceiver);
}

void
MPICommunication::send(int itemToSend, int rankReceiver) {
  preciceTrace2("send(int)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion(rankReceiver != ANY_SENDER);
  MPI_Send(&itemToSend,
           1,
           MPI_INT,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

Request::SharedPointer
MPICommunication::aSend(int* itemToSend, int rankReceiver) {
  return aSend(itemToSend, 1, rankReceiver);
}

void
MPICommunication::send(bool itemToSend, int rankReceiver) {
  preciceTrace2("send(bool)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion(rankReceiver != ANY_SENDER);
  MPI_Send(&itemToSend,
           1,
           MPI_BOOL,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

Request::SharedPointer
MPICommunication::aSend(bool* itemToSend, int rankReceiver) {
  preciceTrace("aSend(bool*)");
  rankReceiver = rankReceiver - _rankOffset;
  assertion(rankReceiver != ANY_SENDER);

  MPI_Request request;

  MPI_Isend(itemToSend,
            1,
            MPI_BOOL,
            rank(rankReceiver),
            0,
            communicator(rankReceiver),
            &request);

  return Request::SharedPointer(new MPIRequest(request));
}

int
MPICommunication::receive(std::string& itemToReceive, int rankSender) {
  preciceTrace2("receive(string)", itemToReceive, rankSender);
  rankSender = rankSender - _rankOffset;
  int length;
  MPI_Status status;
  rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
  MPI_Probe(rank(rankSender), 0, communicator(rankSender), &status);
  MPI_Get_count(&status, MPI_CHAR, &length);
  rankSender = status.MPI_SOURCE;
  preciceDebug("Stringlength = " << length);
  char cstr[length];
  MPI_Recv(cstr,
           length,
           MPI_CHAR,
           rank(rankSender),
           0,
           communicator(rankSender),
           MPI_STATUS_IGNORE);
  itemToReceive = std::string(cstr);
  preciceDebug("Received \"" << itemToReceive << "\" from rank " << rankSender);
  return rankSender;
}

int
MPICommunication::receive(int* itemsToReceive, int size, int rankSender) {
  preciceTrace1("receive(int*)", size);
  rankSender = rankSender - _rankOffset;
  rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
  MPI_Status status;
  MPI_Recv(itemsToReceive,
           size,
           MPI_INT,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
  return status.MPI_SOURCE;
}

int
MPICommunication::receive(double* itemsToReceive, int size, int rankSender) {
  preciceTrace1("receive(double*)", size);
  rankSender = rankSender - _rankOffset;
  rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
  MPI_Status status;
  MPI_Recv(itemsToReceive,
           size,
           MPI_DOUBLE,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
  return status.MPI_SOURCE;
}

int
MPICommunication::receive(double& itemToReceive, int rankSender) {
  preciceTrace1("receive(double)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
  MPI_Status status;
  MPI_Recv(&itemToReceive,
           1,
           MPI_DOUBLE,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
  preciceDebug("Received " << itemToReceive << " from rank "
                           << status.MPI_SOURCE);
  return status.MPI_SOURCE;
}

int
MPICommunication::receive(int& itemToReceive, int rankSender) {
  preciceTrace1("receive(int)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
  MPI_Status status;
  MPI_Recv(&itemToReceive,
           1,
           MPI_INT,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
  preciceDebug("Received " << itemToReceive << " from rank "
                           << status.MPI_SOURCE);
  return status.MPI_SOURCE;
}

int
MPICommunication::receive(bool& itemToReceive, int rankSender) {
  preciceTrace1("receive(bool)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = rankSender == ANY_SENDER ? MPI_ANY_SOURCE : rankSender;
  MPI_Status status;
  MPI_Recv(&itemToReceive,
           1,
           MPI_BOOL,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
  preciceDebug("Received " << itemToReceive << " from rank "
                           << status.MPI_SOURCE);
  return status.MPI_SOURCE;
}
}
} // namespace precice, com

#endif // not PRECICE_NO_MPI
