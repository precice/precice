// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#include "MPIDirectCommunication.hpp"

#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"

#include <map>

namespace precice {
namespace com {

tarch::logging::Log MPIDirectCommunication::_log(
    "precice::com::MPIDirectCommunication");

MPIDirectCommunication::MPIDirectCommunication()
    : _communicator(utils::Parallel::getGlobalCommunicator())
    , _globalCommunicator(utils::Parallel::getGlobalCommunicator())
    , _localCommunicator(utils::Parallel::getLocalCommunicator())
    , _isConnected(false) {
}

MPIDirectCommunication::~MPIDirectCommunication() {
  preciceTrace1("~MPIDirectCommunication()", _isConnected);

  closeConnection();
}

int
MPIDirectCommunication::getRemoteCommunicatorSize() {
  preciceTrace("getRemoteCommunicatorSize()");
  assertion(isConnected());
  int remoteSize = 0;
  MPI_Comm_remote_size(communicator(), &remoteSize);
  return remoteSize;
}

void
MPIDirectCommunication::acceptConnection(std::string const& nameAcceptor,
                                         std::string const& nameRequester,
                                         int acceptorProcessRank,
                                         int acceptorCommunicatorSize) {
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);
  assertion(not isConnected());

  int argc = 1;
  char* arg = new char[8];
  strcpy(arg, "precice");
  char** argv = &arg;
  utils::Parallel::initialize(&argc, &argv, nameAcceptor);
  delete[] arg;

  preciceCheck(utils::Parallel::getCommunicatorSize() > 1,
               "acceptConnection()",
               "ERROR: MPI communication direct (i.e. single) can be only "
                   << "used with more than one process in base communicator!");
  _globalCommunicator = utils::Parallel::getGlobalCommunicator();
  _localCommunicator = utils::Parallel::getLocalCommunicator();
  MPI_Intercomm_create(
      _localCommunicator,
      0, // Local communicator, local leader rank
      _globalCommunicator,
      getLeaderRank(nameRequester), // Peer communicator, remote leader rank
      0,
      &communicator()); // Tag, intercommunicator to be created
  _isConnected = true;
}

void
MPIDirectCommunication::closeConnection() {
  preciceTrace("closeConnection()");

  if (not isConnected())
    return;

  MPI_Comm_free(&communicator());

  _isConnected = false;
}

void
MPIDirectCommunication::requestConnection(std::string const& nameAcceptor,
                                          std::string const& nameRequester,
                                          int requesterProcessRank,
                                          int requesterCommunicatorSize) {
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);
  assertion(not isConnected());

  int argc = 1;
  char* arg = new char[8];
  strcpy(arg, "precice");
  char** argv = &arg;
  utils::Parallel::initialize(&argc, &argv, nameRequester);
  delete[] arg;

  preciceCheck(utils::Parallel::getCommunicatorSize() > 1,
               "requestConnection()",
               "ERROR: MPI communication direct (i.e. single) can be only "
                   << "used with more than one process in base communicator!");
  _globalCommunicator = utils::Parallel::getGlobalCommunicator();
  _localCommunicator = utils::Parallel::getLocalCommunicator();
  MPI_Intercomm_create(
      _localCommunicator,
      0, // Local communicator, local leader rank
      _globalCommunicator,
      getLeaderRank(nameAcceptor), // Peer communicator, remote leader rank
      0,
      &communicator()); // Tag, intercommunicator to be created
  _isConnected = true;
}

int
MPIDirectCommunication::getGroupID(std::string const& accessorName) {
  preciceTrace1("getGroupID()", accessorName);
  typedef utils::Parallel Par;
  const std::vector<Par::AccessorGroup>& _groups = Par::getAccessorGroups();
  for (const Par::AccessorGroup& group : _groups) {
    if (group.name == accessorName) {
      preciceDebug("return group ID = " << group.id);
      return group.id;
    }
  }
  preciceError("getGroupID()",
               "Unknown accessor name \"" << accessorName << "\"!");
}

int
MPIDirectCommunication::getLeaderRank(std::string const& accessorName) {
  preciceTrace1("getLeaderRank()", accessorName);
  typedef utils::Parallel Par;
  const std::vector<Par::AccessorGroup>& _groups = Par::getAccessorGroups();
  for (const Par::AccessorGroup& group : _groups) {
    if (group.name == accessorName) {
      preciceDebug("return rank = " << group.leaderRank);
      return group.leaderRank;
    }
  }
  preciceError("getLeaderRank()",
               "Unknown accessor name \"" << accessorName << "\"!");
}

MPI_Comm&
MPIDirectCommunication::communicator(int rank) {
  return _communicator;
}

int
MPIDirectCommunication::rank(int rank) {
  return rank;
}
}
} // close namespaces

#endif // not PRECICE_NO_MPI
