// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NO_MPI

#include "MPIPortsCommunication.hpp"

#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"

#include <fstream>

namespace precice {
namespace com {
tarch::logging::Log MPIPortsCommunication::_log(
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
  preciceTrace1("~MPIPortsCommunication()", _isConnected);

  if (_isConnected) {
    closeConnection();
  }
}

bool
MPIPortsCommunication::isConnected() {
  return _isConnected;
}

int
MPIPortsCommunication::getRemoteCommunicatorSize() {
  preciceTrace("getRemoteCommunicatorSize()");

  assertion(_isConnected);

  return _communicators.size();
}

void
MPIPortsCommunication::acceptConnection(std::string const& nameAcceptor,
                                        std::string const& nameRequester,
                                        int acceptorProcessRank,
                                        int acceptorCommunicatorSize) {
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);
  preciceCheck(acceptorCommunicatorSize == 1,
               "acceptConnection()",
               "Acceptor of MPI port connection can only have one process!");

  assertion(not _isConnected);

  _isAcceptor = true;

  // BUG:
  // It is extremely important that the call to `Parallel::initialize' follows
  // *after* the call to `MPI_Open_port'. Otherwise, on Windows, even with the
  // latest Intel MPI, the program hangs. Possibly `Parallel::initialize' is
  // doing something weird inside?

  MPI_Open_port(MPI_INFO_NULL, _portName);

  utils::Parallel::initialize(NULL, NULL, nameAcceptor);

  std::string addressFileName(_addressDirectory + "/" + "." + nameRequester +
                              "-" + nameAcceptor + ".address");

  {
    std::ofstream addressFile(addressFileName + "~", std::ios::out);

    addressFile << _portName;
  }

  std::rename((addressFileName + "~").c_str(), addressFileName.c_str());

  MPI_Status status;
  MPI_Comm communicator;

  int requesterProcessRank = -1;
  int requesterCommunicatorSize = 0;

  MPI_Comm_accept(_portName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);

  MPI_Recv(&requesterProcessRank, 1, MPI_INT, 0, 42, communicator, &status);
  MPI_Recv(
      &requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator, &status);

  preciceCheck(requesterCommunicatorSize > 0,
               "acceptConnection()",
               "Requester communicator size has to be > 0!");

  _communicators.resize(requesterCommunicatorSize, MPI_COMM_NULL);

  _communicators[requesterProcessRank] = communicator;

  for (int i = 1; i < requesterCommunicatorSize; ++i) {
    MPI_Comm_accept(_portName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);

    MPI_Recv(&requesterProcessRank, 1, MPI_INT, 0, 42, communicator, &status);
    MPI_Recv(
        &requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator, &status);

    preciceCheck(requesterCommunicatorSize == _communicators.size(),
                 "acceptConnection()",
                 "Requester communicator sizes are inconsistent!");

    preciceCheck(_communicators[requesterProcessRank] == MPI_COMM_NULL,
                 "acceptConnection()",
                 "Duplicate request to connect by same rank ("
                     << requesterProcessRank << ")!");

    _communicators[requesterProcessRank] = communicator;
  }

  std::remove(addressFileName.c_str());

  _isConnected = true;
}

void
MPIPortsCommunication::requestConnection(std::string const& nameAcceptor,
                                         std::string const& nameRequester,
                                         int requesterProcessRank,
                                         int requesterCommunicatorSize) {
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);

  assertion(not _isConnected);

  _isAcceptor = false;

  utils::Parallel::initialize(NULL, NULL, nameRequester);

  std::string addressFileName(_addressDirectory + "/" + "." + nameRequester +
                              "-" + nameAcceptor + ".address");

  {
    std::ifstream addressFile;

    do {
      addressFile.open(addressFileName, std::ios::in);
    } while (not addressFile);

    addressFile.getline(_portName, MPI_MAX_PORT_NAME);
  }

  // std::cout << _portName << std::endl;

  MPI_Comm communicator;

  MPI_Comm_connect(_portName, MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);

  MPI_Send(&requesterProcessRank, 1, MPI_INT, 0, 42, communicator);
  MPI_Send(&requesterCommunicatorSize, 1, MPI_INT, 0, 42, communicator);

  _communicators.resize(1, MPI_COMM_NULL);

  _communicators[0] = communicator;

  _isConnected = true;
}

void
MPIPortsCommunication::closeConnection() {
  assertion(_isConnected);

  for (auto communicator : _communicators) {
    MPI_Comm_disconnect(&communicator);
  }

  if (_isAcceptor) {
    MPI_Close_port(_portName);
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
