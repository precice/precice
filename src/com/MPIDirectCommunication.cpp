#include <cstddef>
#ifndef PRECICE_NO_MPI

#include <memory>

#include "MPIDirectCommunication.hpp"
#include "logging/LogMacros.hpp"
#include "precice/types.hpp"
#include "utils/Parallel.hpp"
#include "utils/assertion.hpp"
#include "utils/span_tools.hpp"

namespace precice {
namespace com {
MPIDirectCommunication::MPIDirectCommunication()
    : _commState(utils::Parallel::current())
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
                                              std::string const &tag,
                                              int                acceptorRank,
                                              int                rankOffset)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected());
  // MPI Direct Comm only supports MasterSlave connections
  PRECICE_ASSERT(rankOffset == 1, "MPIDirectCommunication only supports MasterSlave Communications!");
  setRankOffset(rankOffset);

  _commState   = utils::Parallel::current();
  _isConnected = true;

  PRECICE_ASSERT(acceptorRank == 0, "The Acceptor/Master has to be rank 0!");
  PRECICE_ASSERT(_commState->rank() == acceptorRank, "The given acceptor rank does not match the communicator rank!");
}

void MPIDirectCommunication::closeConnection()
{
  PRECICE_TRACE();

  if (not isConnected())
    return;

  _isConnected = false;
}

void MPIDirectCommunication::requestConnection(std::string const &acceptorName,
                                               std::string const &requesterName,
                                               std::string const &tag,
                                               int                requesterRank,
                                               int                requesterCommunicatorSize)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected());

  setRankOffset(0); //rankOffset makes no sense here
  _commState   = utils::Parallel::current();
  _isConnected = true;

  PRECICE_ASSERT(requesterRank == _commState->rank() - 1);
  PRECICE_ASSERT(requesterCommunicatorSize + 1 == _commState->size());
}

void MPIDirectCommunication::reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());
  Rank rank = _commState->rank();
  MPI_Reduce(const_cast<double *>(itemsToSend.data()), itemsToReceive.data(), itemsToSend.size(), MPI_DOUBLE, MPI_SUM, rank, _commState->comm);
}

void MPIDirectCommunication::reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank rankMaster)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());
  MPI_Reduce(const_cast<double *>(itemsToSend.data()), itemsToReceive.data(), itemsToSend.size(), MPI_DOUBLE, MPI_SUM, rankMaster, _commState->comm);
}

void MPIDirectCommunication::reduceSum(int itemToSend, int &itemsToReceive)
{
  PRECICE_TRACE();
  Rank rank = _commState->rank();
  MPI_Reduce(&itemToSend, &itemsToReceive, 1, MPI_INT, MPI_SUM, rank, _commState->comm);
}

void MPIDirectCommunication::reduceSum(int itemToSend, int &itemsToReceive, Rank rankMaster)
{
  PRECICE_TRACE();
  MPI_Reduce(&itemToSend, &itemsToReceive, 1, MPI_INT, MPI_SUM, rankMaster, _commState->comm);
}

void MPIDirectCommunication::allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());
  MPI_Allreduce(const_cast<double *>(itemsToSend.data()), itemsToReceive.data(), itemsToSend.size(), MPI_DOUBLE, MPI_SUM, _commState->comm);
}

void MPIDirectCommunication::allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank rankMaster)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());
  MPI_Allreduce(const_cast<double *>(itemsToSend.data()), itemsToReceive.data(), itemsToReceive.size(), MPI_DOUBLE, MPI_SUM, _commState->comm);
}

void MPIDirectCommunication::allreduceSum(double itemToSend, double &itemToReceive)
{
  PRECICE_TRACE();
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_DOUBLE, MPI_SUM, _commState->comm);
}

void MPIDirectCommunication::allreduceSum(double itemToSend, double &itemToReceive, Rank rankMaster)
{
  PRECICE_TRACE();
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_DOUBLE, MPI_SUM, _commState->comm);
}

void MPIDirectCommunication::allreduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_INT, MPI_SUM, _commState->comm);
}

void MPIDirectCommunication::allreduceSum(int itemToSend, int &itemToReceive, Rank rankMaster)
{
  PRECICE_TRACE();
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_INT, MPI_SUM, _commState->comm);
}

void MPIDirectCommunication::broadcast(precice::span<const int> itemsToSend)
{
  PRECICE_TRACE(itemsToSend.size());
  MPI_Bcast(const_cast<int *>(itemsToSend.data()), itemsToSend.size(), MPI_INT, 0, _commState->comm);
}

void MPIDirectCommunication::broadcast(precice::span<int> itemsToReceive, int rankBroadcaster)
{
  PRECICE_TRACE(itemsToReceive.size());
  MPI_Bcast(itemsToReceive.data(), itemsToReceive.size(), MPI_INT, rankBroadcaster, _commState->comm);
}

void MPIDirectCommunication::broadcast(int itemToSend)
{
  PRECICE_TRACE();
  broadcast(precice::refToSpan<const int>(itemToSend));
}

void MPIDirectCommunication::broadcast(int &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  broadcast(precice::refToSpan<int>(itemToReceive), rankBroadcaster);
}

void MPIDirectCommunication::broadcast(precice::span<const double> itemsToSend)
{
  PRECICE_TRACE(itemsToSend.size());
  MPI_Bcast(const_cast<double *>(itemsToSend.data()), itemsToSend.size(), MPI_DOUBLE, 0, _commState->comm);
}

void MPIDirectCommunication::broadcast(precice::span<double> itemsToReceive, int rankBroadcaster)
{
  PRECICE_TRACE(itemsToReceive.size());
  MPI_Bcast(itemsToReceive.data(), itemsToReceive.size(), MPI_DOUBLE, rankBroadcaster, _commState->comm);
}

void MPIDirectCommunication::broadcast(double itemToSend)
{
  PRECICE_TRACE();
  broadcast(precice::refToSpan<const double>(itemToSend));
}

void MPIDirectCommunication::broadcast(double &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  broadcast(precice::refToSpan<double>(itemToReceive), rankBroadcaster);
}

void MPIDirectCommunication::broadcast(bool itemToSend)
{
  PRECICE_TRACE();
  int item = itemToSend;
  broadcast(item);
}

void MPIDirectCommunication::broadcast(bool &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  int item;
  broadcast(item, rankBroadcaster);
  itemToReceive = item;
}

MPI_Comm &MPIDirectCommunication::communicator(Rank rank)
{
  return _commState->comm;
}

int MPIDirectCommunication::rank(Rank rank)
{
  // Correct _rankOffset if we are on master
  return rank;
}

int MPIDirectCommunication::adjustRank(Rank rank) const
{
  return rank;
}

} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
