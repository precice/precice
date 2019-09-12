#ifndef PRECICE_NO_MPI

#include "MPIDirectCommunication.hpp"
#include "utils/Parallel.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace com {
MPIDirectCommunication::MPIDirectCommunication()
    : _communicator(utils::Parallel::getGlobalCommunicator()),
      _globalCommunicator(utils::Parallel::getGlobalCommunicator()),
      _localCommunicator(utils::Parallel::getGlobalCommunicator())
{
}

MPIDirectCommunication::~MPIDirectCommunication()
{
  PRECICE_TRACE(_isConnected);
  closeConnection();
}

size_t MPIDirectCommunication::getRemoteCommunicatorSize()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(isConnected());
  int remoteSize = 0;
  MPI_Comm_remote_size(communicator(), &remoteSize);
  return remoteSize;
}

void MPIDirectCommunication::acceptConnection(std::string const &acceptorName,
                                              std::string const &requesterName,
                                              int                acceptorRank)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected());

  utils::Parallel::splitCommunicator(acceptorName);

  PRECICE_CHECK(utils::Parallel::getCommunicatorSize() > 1,
                "MPI communication direct (i.e. single) can be only used with more than one process in base communicator!");

  _globalCommunicator = utils::Parallel::getGlobalCommunicator();
  _localCommunicator  = utils::Parallel::getLocalCommunicator();
  MPI_Intercomm_create(
      _localCommunicator,
      0, // Local communicator, local leader rank
      _globalCommunicator,
      getLeaderRank(requesterName), // Peer communicator, remote leader rank
      0,
      &communicator()); // Tag, intercommunicator to be created
  _isConnected = true;
}

void MPIDirectCommunication::closeConnection()
{
  PRECICE_TRACE()

  if (not isConnected())
    return;

  MPI_Comm_free(&communicator());
  _isConnected = false;
}

void MPIDirectCommunication::requestConnection(std::string const &acceptorName,
                                               std::string const &requesterName,
                                               int                requesterRank,
                                               int                requesterCommunicatorSize)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected());

  utils::Parallel::splitCommunicator(requesterName);

  PRECICE_CHECK(utils::Parallel::getCommunicatorSize() > 1,
                "MPI communication direct (i.e. single) can be only used with more than one process in base communicator!");

  _globalCommunicator = utils::Parallel::getGlobalCommunicator();
  _localCommunicator  = utils::Parallel::getLocalCommunicator();
  MPI_Intercomm_create(
      _localCommunicator,
      0, // Local communicator, local leader rank
      _globalCommunicator,
      getLeaderRank(acceptorName), // Peer communicator, remote leader rank
      0,
      &communicator()); // Tag, intercommunicator to be created
  _isConnected = true;
}

int MPIDirectCommunication::getGroupID(std::string const &accessorName)
{
  PRECICE_TRACE(accessorName);
  using Par                                      = utils::Parallel;
  const std::vector<Par::AccessorGroup> &_groups = Par::getAccessorGroups();
  for (const Par::AccessorGroup &group : _groups) {
    if (group.name == accessorName) {
      PRECICE_DEBUG("return group ID = " << group.id);
      return group.id;
    }
  }
  PRECICE_ERROR("Unknown accessor name \"" << accessorName << "\"!");
}

int MPIDirectCommunication::getLeaderRank(std::string const &accessorName)
{
  PRECICE_TRACE(accessorName);
  using Par                                      = utils::Parallel;
  const std::vector<Par::AccessorGroup> &_groups = Par::getAccessorGroups();
  for (const Par::AccessorGroup &group : _groups) {
    if (group.name == accessorName) {
      PRECICE_DEBUG("return rank = " << group.leaderRank);
      return group.leaderRank;
    }
  }
  PRECICE_ERROR("Unknown accessor name \"" << accessorName << "\"!");
}

void MPIDirectCommunication::reduceSum(double *itemsToSend, double *itemsToReceive, int size)
{
  PRECICE_TRACE(size);
  int rank = -1;
  MPI_Comm_rank(_globalCommunicator, &rank);
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Reduce(itemsToSend, itemsToReceive, size, MPI_DOUBLE, MPI_SUM, rank, _globalCommunicator);
}

void MPIDirectCommunication::reduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster)
{
  PRECICE_TRACE(size);
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Reduce(itemsToSend, itemsToReceive, size, MPI_DOUBLE, MPI_SUM, rankMaster, _globalCommunicator);
}

void MPIDirectCommunication::reduceSum(int itemToSend, int &itemsToReceive)
{
  PRECICE_TRACE();
  int rank = -1;
  MPI_Comm_rank(_globalCommunicator, &rank);
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Reduce(&itemToSend, &itemsToReceive, 1, MPI_INT, MPI_SUM, rank, _globalCommunicator);
}

void MPIDirectCommunication::reduceSum(int itemToSend, int &itemsToReceive, int rankMaster)
{
  PRECICE_TRACE();
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Reduce(&itemToSend, &itemsToReceive, 1, MPI_INT, MPI_SUM, rankMaster, _globalCommunicator);
}

void MPIDirectCommunication::allreduceSum(double *itemsToSend, double *itemsToReceive, int size)
{
  PRECICE_TRACE(size);
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Allreduce(itemsToSend, itemsToReceive, size, MPI_DOUBLE, MPI_SUM, _globalCommunicator);
}

void MPIDirectCommunication::allreduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster)
{
  PRECICE_TRACE(size);
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Allreduce(itemsToSend, itemsToReceive, size, MPI_DOUBLE, MPI_SUM, _globalCommunicator);
}

void MPIDirectCommunication::allreduceSum(double itemToSend, double &itemToReceive)
{
  PRECICE_TRACE();
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_DOUBLE, MPI_SUM, _globalCommunicator);
}

void MPIDirectCommunication::allreduceSum(double itemToSend, double &itemToReceive, int rankMaster)
{
  PRECICE_TRACE();
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_DOUBLE, MPI_SUM, _globalCommunicator);
}

void MPIDirectCommunication::allreduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_INT, MPI_SUM, _globalCommunicator);
}

void MPIDirectCommunication::allreduceSum(int itemToSend, int &itemToReceive, int rankMaster)
{
  PRECICE_TRACE();
  // _comunicator did't work here as we seem to have two communicators, one with the master and one with the slaves
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_INT, MPI_SUM, _globalCommunicator);
}

void MPIDirectCommunication::broadcast(const int *itemsToSend, int size)
{
  PRECICE_TRACE(size);
  MPI_Bcast(const_cast<int *>(itemsToSend), size, MPI_INT, MPI_ROOT, _communicator);
}

void MPIDirectCommunication::broadcast(int *itemsToReceive,
                                       int  size,
                                       int  rankBroadcaster)
{
  PRECICE_TRACE(size);
  MPI_Bcast(itemsToReceive, size, MPI_INT, rankBroadcaster, _communicator);
}

void MPIDirectCommunication::broadcast(int itemToSend)
{
  PRECICE_TRACE();
  broadcast(&itemToSend, 1);
}

void MPIDirectCommunication::broadcast(int &itemToReceive, int rankBroadcaster)
{
  PRECICE_TRACE();
  broadcast(&itemToReceive, 1, rankBroadcaster);
}

void MPIDirectCommunication::broadcast(const double *itemsToSend, int size)
{
  PRECICE_TRACE(size);
  MPI_Bcast(const_cast<double *>(itemsToSend), size, MPI_DOUBLE, MPI_ROOT, _communicator);
}

void MPIDirectCommunication::broadcast(double *itemsToReceive,
                                       int     size,
                                       int     rankBroadcaster)
{
  PRECICE_TRACE(size);
  MPI_Bcast(itemsToReceive, size, MPI_DOUBLE, rankBroadcaster, _communicator);
}

void MPIDirectCommunication::broadcast(double itemToSend)
{
  PRECICE_TRACE();
  broadcast(&itemToSend, 1);
}

void MPIDirectCommunication::broadcast(double &itemToReceive, int rankBroadcaster)
{
  PRECICE_TRACE();
  broadcast(&itemToReceive, 1, rankBroadcaster);
}

void MPIDirectCommunication::broadcast(bool itemToSend)
{
  PRECICE_TRACE();
  int item = itemToSend;
  broadcast(item);
}

void MPIDirectCommunication::broadcast(bool &itemToReceive, int rankBroadcaster)
{
  PRECICE_TRACE();
  int item;
  broadcast(item, rankBroadcaster);
  itemToReceive = item;
}

MPI_Comm &MPIDirectCommunication::communicator(int rank)
{
  return _communicator;
}

int MPIDirectCommunication::rank(int rank)
{
  return rank;
}
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
