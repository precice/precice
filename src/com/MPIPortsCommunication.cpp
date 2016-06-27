// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NO_MPI

#include "MPIPortsCommunication.hpp"

#include "utils/Globals.hpp"
#include "utils/Publisher.hpp"
#include "utils/Parallel.hpp"

#include <chrono>
#include <sstream>
#include <thread>

using precice::utils::Publisher;
using precice::utils::ScopedPublisher;

namespace precice {
namespace com {
logging::Logger MPIPortsCommunication::_log(
    "precice::com::MPIPortsCommunication");

MPIPortsCommunication::MPIPortsCommunication(
    std::string const& addressDirectory)
    : _addressDirectory(addressDirectory)
    , _isAcceptor(false)
    , _isConnected(false) {
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

MPIPortsCommunication::~MPIPortsCommunication() {
  preciceTrace("~MPIPortsCommunication()", _isConnected);

  closeConnection();
}

bool
MPIPortsCommunication::isConnected() {
  return _isConnected;
}

size_t MPIPortsCommunication::getRemoteCommunicatorSize() {
  preciceTrace("getRemoteCommunicatorSize()");

  assertion(isConnected());

  return _communicators.size();
}

void
MPIPortsCommunication::acceptConnection(std::string const& nameAcceptor,
                                        std::string const& nameRequester,
                                        int acceptorProcessRank,
                                        int acceptorCommunicatorSize) {
  preciceTrace("acceptConnection()", nameAcceptor, nameRequester);

  preciceCheck(acceptorCommunicatorSize == 1,
               "acceptConnection()",
               "Acceptor of MPI port connection can only have one process!");

  assertion(not isConnected());

  _isAcceptor = true;
  _rank = acceptorProcessRank;

  // BUG report from Alex:
  // It is extremely important that the call to `Parallel::initialize' follows
  // *after* the call to `MPI_Open_port'. Otherwise, on Windows, even with the
  // latest Intel MPI, the program hangs. Possibly `Parallel::initialize' is
  // doing something weird inside?

  MPI_Open_port(MPI_INFO_NULL, _portName);

  std::string address(_portName);
  std::string addressFileName("." + nameRequester + "-" + nameAcceptor +
                              ".address");

  Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);

  ScopedPublisher p(addressFileName);

  p.write(address);

  preciceDebug("Accept connection at " << address);

  MPI_Comm communicator;

  MPI_Comm_accept(_portName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);

  preciceDebug("Accepted connection at " << address);

  int requesterProcessRank = -1;
  int requesterCommunicatorSize = 0;

  MPI_Recv(&requesterProcessRank,
           1,
           MPI_INT,
           0,
           42,
           communicator,
           MPI_STATUS_IGNORE);
  MPI_Recv(&requesterCommunicatorSize,
           1,
           MPI_INT,
           0,
           42,
           communicator,
           MPI_STATUS_IGNORE);

  preciceCheck(requesterCommunicatorSize > 0,
               "acceptConnection()",
               "Requester communicator "
                   << "size has to be > 0!");

  _communicators.resize(requesterCommunicatorSize, MPI_COMM_NULL);

  _communicators[requesterProcessRank] = communicator;

  _isConnected = true;

  for (int i = 1; i < requesterCommunicatorSize; ++i) {
    MPI_Comm_accept(_portName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);

    preciceDebug("Accepted connection at " << address);

    MPI_Recv(&requesterProcessRank,
             1,
             MPI_INT,
             0,
             42,
             communicator,
             MPI_STATUS_IGNORE);
    MPI_Recv(&requesterCommunicatorSize,
             1,
             MPI_INT,
             0,
             42,
             communicator,
             MPI_STATUS_IGNORE);

    preciceCheck(requesterCommunicatorSize == _communicators.size(),
                 "acceptConnection()",
                 "Requester communicator sizes are inconsistent!");
    preciceCheck(_communicators[requesterProcessRank] == MPI_COMM_NULL,
                 "acceptConnection()",
                 "Duplicate request to connect by same rank ("
                     << requesterProcessRank << ")!");

    _communicators[requesterProcessRank] = communicator;

    _isConnected = true;
  }
}

void
MPIPortsCommunication::acceptConnectionAsServer(
    std::string const& nameAcceptor,
    std::string const& nameRequester,
    int requesterCommunicatorSize) {
  preciceTrace("acceptConnectionAsServer()", nameAcceptor, nameRequester);

  preciceCheck(requesterCommunicatorSize > 0,
               "acceptConnectionAsServer()",
               "Requester communicator "
                   << "size has to be > 0!");

  assertion(not isConnected());

  _isAcceptor = true;
  _rank = 0;

  // BUG report from Alex:
  // It is extremely important that the call to `Parallel::initialize' follows
  // *after* the call to `MPI_Open_port'. Otherwise, on Windows, even with the
  // latest Intel MPI, the program hangs. Possibly `Parallel::initialize' is
  // doing something weird inside?

  MPI_Open_port(MPI_INFO_NULL, _portName);

  std::string address(_portName);
  std::string addressFileName("." + nameRequester + "-" + nameAcceptor +
                              ".address");

  Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);

  ScopedPublisher p(addressFileName);

  p.write(address);

  preciceDebug("Accept connection at " << address);

  _communicators.resize(requesterCommunicatorSize, MPI_COMM_NULL);

  for (int requesterProcessRank = 0;
       requesterProcessRank < requesterCommunicatorSize;
       ++requesterProcessRank) {
    MPI_Comm communicator;

    MPI_Comm_accept(_portName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);

    preciceDebug("Accepted connection at " << address);

    preciceCheck(_communicators[requesterProcessRank] == MPI_COMM_NULL,
                 "acceptConnectionAsServer()",
                 "Duplicate request to connect by same rank ("
                     << requesterProcessRank << ")!");

    _communicators[requesterProcessRank] = communicator;

    _isConnected = true;

    // BUG:
    // On Windows, with Intel MPI, in point-to-point communication integration
    // test, a deadlock happens because the following `MPI_Send' has no effect.
    // This happens rarely and the nature of this phenomenon is still
    // unknown. It looks as if the actual message (to be sent) is lost
    // somehow. To me that one looks more like an implementation bug. Let's see
    // how it goes in other environments (OS/MPI combinations).
    MPI_Send(&requesterProcessRank, 1, MPI_INT, 0, 42, communicator);
    // send(requesterProcessRank, requesterProcessRank);
  }
}

void
MPIPortsCommunication::requestConnection(std::string const& nameAcceptor,
                                         std::string const& nameRequester,
                                         int requesterProcessRank,
                                         int requesterCommunicatorSize) {
  preciceTrace("requestConnection()", nameAcceptor, nameRequester);

  assertion(not isConnected());

  _isAcceptor = false;

  std::string address;
  std::string addressFileName("." + nameRequester + "-" + nameAcceptor +
                              ".address");

  Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);

  Publisher p(addressFileName);

  p.read(address);

  preciceDebug("Request connection to " << address);

  {
    std::istringstream iss(address);

    iss >> _portName;
  }

  MPI_Comm communicator;

  MPI_Comm_connect(_portName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);

  preciceDebug("Requested connection to " << address);

  _communicators.push_back(communicator);

  _isConnected = true;

  _rank = requesterProcessRank;

  MPI_Send(&requesterProcessRank, 1, MPI_INT, 0, 42, communicator);
  MPI_Send(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator);
}

int
MPIPortsCommunication::requestConnectionAsClient(
    std::string const& nameAcceptor, std::string const& nameRequester) {
  preciceTrace("requestConnectionAsClient()", nameAcceptor, nameRequester);

  assertion(not isConnected());

  _isAcceptor = false;

  std::string address;
  std::string addressFileName("." + nameRequester + "-" + nameAcceptor +
                              ".address");

  Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);

  Publisher p(addressFileName);

  p.read(address);

  preciceDebug("Request connection to " << address);

  {
    std::istringstream iss(address);

    iss >> _portName;
  }

  MPI_Comm communicator;

  MPI_Comm_connect(_portName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);

  preciceDebug("Requested connection to " << address);

  _communicators.push_back(communicator);

  _isConnected = true;

  // BUG:
  // On Windows, with Intel MPI, in point-to-point communication integration
  // test, a deadlock happens because the following `MPI_Recv' never
  // returns. This happens rarely and the nature of this phenomenon is still
  // unknown. It looks as if the actual message (to be received) was lost
  // somehow. To me that one looks more like an implementation bug. Let's see
  // how it goes in other environments (OS/MPI combinations).
  // MPI_Recv(&_rank, 1, MPI_INT, 0, 42, communicator, MPI_STATUS_IGNORE);
  // receive(_rank, 0);

  // NOTE:
  // This is a partial solution. First of all, somehow it seems to lower the
  // deadlock frequency. Secondly, it gives some time to ensure that there is a
  // deadlock. Finally, when the deadlock really happens, it reports a proper
  // error and terminates the application.
  {
    MPI_Request request;

    MPI_Irecv(&_rank, 1, MPI_INT, 0, 42, communicator, &request);

    int complete = 0;
    int i;

    for (i = 0; not complete && i < 500; ++i) {
      MPI_Test(&request, &complete, MPI_STATUS_IGNORE);

      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    if (i >= 500) {
      preciceError("requestConnectionAsClient()",
                   "Oops, we have a deadlock here..."
                   " "
                   "Now terminating, retry please!");

      exit(1);
    }
  }

  return _rank;
}

void
MPIPortsCommunication::closeConnection() {
  preciceTrace("closeConnection()", _communicators.size());

  if (not isConnected())
    return;

  for (auto communicator : _communicators) {
    MPI_Comm_disconnect(&communicator);
  }

  preciceDebug("Disconnected");

  if (_isAcceptor) {
    MPI_Close_port(_portName);
    preciceDebug("Port closed");
  }

  _isConnected = false;
}

MPI_Comm&
MPIPortsCommunication::communicator(int rank) {
  return _communicators[rank];
}

int
MPIPortsCommunication::rank(int rank) {
  return 0;
}
}
} // namespace precice, com

#endif // not PRECICE_NO_MPI
