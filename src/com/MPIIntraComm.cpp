#ifndef PRECICE_NO_MPI

#include "com/MPIIntraComm.hpp"
#include "com/MPIRequest.hpp"
#include "logging/LogMacros.hpp"
#include "precice/impl/Types.hpp"
#include "utils/Parallel.hpp"
#include "utils/assertion.hpp"
#include "utils/span_tools.hpp"

#include <cstddef>
#include <memory>
#include <ostream>

template <size_t>
struct MPI_Select_unsigned_integer_datatype_intra;

template <>
struct MPI_Select_unsigned_integer_datatype_intra<1> {
  static MPI_Datatype datatype;
};
MPI_Datatype MPI_Select_unsigned_integer_datatype_intra<1>::datatype = MPI_UNSIGNED_CHAR;

template <>
struct MPI_Select_unsigned_integer_datatype_intra<2> {
  static MPI_Datatype datatype;
};
MPI_Datatype MPI_Select_unsigned_integer_datatype_intra<2>::datatype = MPI_UNSIGNED_SHORT;

template <>
struct MPI_Select_unsigned_integer_datatype_intra<4> {
  static MPI_Datatype datatype;
};
MPI_Datatype MPI_Select_unsigned_integer_datatype_intra<4>::datatype = MPI_UNSIGNED;

template <>
struct MPI_Select_unsigned_integer_datatype_intra<8> {
  static MPI_Datatype datatype;
};
MPI_Datatype MPI_Select_unsigned_integer_datatype_intra<8>::datatype = MPI_UNSIGNED_LONG;

#define MPI_BOOL_INTRA MPI_Select_unsigned_integer_datatype_intra<sizeof(bool)>::datatype

namespace precice::com {

MPIIntraComm::MPIIntraComm()
    : _commState(utils::Parallel::current())
{
}

MPIIntraComm::~MPIIntraComm()
{
  PRECICE_TRACE(_isConnected);
  closeConnection();
}

bool MPIIntraComm::isConnected()
{
  return _isConnected;
}

size_t MPIIntraComm::getRemoteCommunicatorSize()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(isConnected());
  // In intra-comm, remote = secondaries = size - 1
  return static_cast<size_t>(_commState->size() - 1);
}

void MPIIntraComm::connectIntraComm(std::string const &participantName,
                                    std::string const &tag,
                                    int                rank,
                                    int                size)
{
  PRECICE_TRACE(participantName, rank, size);
  if (size == 1)
    return;

  // For MPI direct, the communicator is already set up by Parallel::splitCommunicator.
  // Just store the current communicator state.
  _commState   = utils::Parallel::current();
  _isConnected = true;

  PRECICE_DEBUG("Connected MPIIntraComm for participant {}, rank {}/{}", participantName, rank, size);
}

void MPIIntraComm::closeConnection()
{
  PRECICE_TRACE();
  if (not isConnected())
    return;
  _isConnected = false;
}

// ==========================================================================
// Reduction — directly using MPI_Reduce / MPI_Allreduce
// ==========================================================================

void MPIIntraComm::reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());
  Rank myRank = _commState->rank();
  MPI_Reduce(const_cast<double *>(itemsToSend.data()), itemsToReceive.data(),
             itemsToSend.size(), MPI_DOUBLE, MPI_SUM, myRank, _commState->comm);
}

void MPIIntraComm::reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());
  MPI_Reduce(const_cast<double *>(itemsToSend.data()), itemsToReceive.data(),
             itemsToSend.size(), MPI_DOUBLE, MPI_SUM, primaryRank, _commState->comm);
}

void MPIIntraComm::reduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();
  Rank myRank = _commState->rank();
  MPI_Reduce(&itemToSend, &itemToReceive, 1, MPI_INT, MPI_SUM, myRank, _commState->comm);
}

void MPIIntraComm::reduceSum(int itemToSend, int &itemToReceive, Rank primaryRank)
{
  PRECICE_TRACE();
  MPI_Reduce(&itemToSend, &itemToReceive, 1, MPI_INT, MPI_SUM, primaryRank, _commState->comm);
}

void MPIIntraComm::allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());
  MPI_Allreduce(const_cast<double *>(itemsToSend.data()), itemsToReceive.data(),
                itemsToSend.size(), MPI_DOUBLE, MPI_SUM, _commState->comm);
}

void MPIIntraComm::allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());
  // MPI_Allreduce is collective, primaryRank parameter is unused but kept for interface compatibility
  MPI_Allreduce(const_cast<double *>(itemsToSend.data()), itemsToReceive.data(),
                itemsToReceive.size(), MPI_DOUBLE, MPI_SUM, _commState->comm);
}

void MPIIntraComm::allreduceSum(double itemToSend, double &itemToReceive)
{
  PRECICE_TRACE();
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_DOUBLE, MPI_SUM, _commState->comm);
}

void MPIIntraComm::allreduceSum(double itemToSend, double &itemToReceive, Rank primaryRank)
{
  PRECICE_TRACE();
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_DOUBLE, MPI_SUM, _commState->comm);
}

void MPIIntraComm::allreduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_INT, MPI_SUM, _commState->comm);
}

void MPIIntraComm::allreduceSum(int itemToSend, int &itemToReceive, Rank primaryRank)
{
  PRECICE_TRACE();
  MPI_Allreduce(&itemToSend, &itemToReceive, 1, MPI_INT, MPI_SUM, _commState->comm);
}

// ==========================================================================
// Broadcast — directly using MPI_Bcast
// ==========================================================================

void MPIIntraComm::broadcast(precice::span<const int> itemsToSend)
{
  PRECICE_TRACE(itemsToSend.size());
  MPI_Bcast(const_cast<int *>(itemsToSend.data()), itemsToSend.size(), MPI_INT, 0, _commState->comm);
}

void MPIIntraComm::broadcast(precice::span<int> itemsToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE(itemsToReceive.size());
  MPI_Bcast(itemsToReceive.data(), itemsToReceive.size(), MPI_INT, rankBroadcaster, _commState->comm);
}

void MPIIntraComm::broadcast(int itemToSend)
{
  PRECICE_TRACE();
  MPI_Bcast(&itemToSend, 1, MPI_INT, 0, _commState->comm);
}

void MPIIntraComm::broadcast(int &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  MPI_Bcast(&itemToReceive, 1, MPI_INT, rankBroadcaster, _commState->comm);
}

void MPIIntraComm::broadcast(precice::span<const double> itemsToSend)
{
  PRECICE_TRACE(itemsToSend.size());
  MPI_Bcast(const_cast<double *>(itemsToSend.data()), itemsToSend.size(), MPI_DOUBLE, 0, _commState->comm);
}

void MPIIntraComm::broadcast(precice::span<double> itemsToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE(itemsToReceive.size());
  MPI_Bcast(itemsToReceive.data(), itemsToReceive.size(), MPI_DOUBLE, rankBroadcaster, _commState->comm);
}

void MPIIntraComm::broadcast(double itemToSend)
{
  PRECICE_TRACE();
  MPI_Bcast(&itemToSend, 1, MPI_DOUBLE, 0, _commState->comm);
}

void MPIIntraComm::broadcast(double &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  MPI_Bcast(&itemToReceive, 1, MPI_DOUBLE, rankBroadcaster, _commState->comm);
}

void MPIIntraComm::broadcast(bool itemToSend)
{
  PRECICE_TRACE();
  int item = itemToSend;
  broadcast(item);
}

void MPIIntraComm::broadcast(bool &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  int item;
  broadcast(item, rankBroadcaster);
  itemToReceive = item;
}

void MPIIntraComm::broadcast(std::vector<int> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(precice::span<const int>{v});
}

void MPIIntraComm::broadcast(std::vector<int> &v, Rank rankBroadcaster)
{
  int size = 0;
  broadcast(size, rankBroadcaster);
  v.clear();
  v.resize(size);
  broadcast(precice::span<int>{v}, rankBroadcaster);
}

void MPIIntraComm::broadcast(std::vector<double> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(precice::span<const double>{v});
}

void MPIIntraComm::broadcast(std::vector<double> &v, Rank rankBroadcaster)
{
  int size = 0;
  broadcast(size, rankBroadcaster);
  v.clear();
  v.resize(size);
  broadcast(precice::span<double>{v}, rankBroadcaster);
}

// ==========================================================================
// Point-to-point — directly using MPI_Send / MPI_Recv / MPI_Isend / MPI_Irecv
// No rank adjustment needed. Ranks are participant ranks directly.
// ==========================================================================

void MPIIntraComm::send(std::string const &itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);
  MPI_Send(const_cast<char *>(itemToSend.c_str()),
           itemToSend.size(), MPI_CHAR, rankReceiver, 0, _commState->comm);
}

void MPIIntraComm::send(precice::span<const int> itemsToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemsToSend.size());
  MPI_Send(const_cast<int *>(itemsToSend.data()),
           itemsToSend.size(), MPI_INT, rankReceiver, 0, _commState->comm);
}

PtrRequest MPIIntraComm::aSend(precice::span<const int> itemsToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemsToSend.size());
  MPI_Request request;
  MPI_Isend(const_cast<int *>(itemsToSend.data()),
            itemsToSend.size(), MPI_INT, rankReceiver, 0, _commState->comm, &request);
  return PtrRequest(new MPIRequest(request));
}

void MPIIntraComm::send(precice::span<const double> itemsToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemsToSend.size());
  MPI_Send(const_cast<double *>(itemsToSend.data()),
           itemsToSend.size(), MPI_DOUBLE, rankReceiver, 0, _commState->comm);
}

PtrRequest MPIIntraComm::aSend(precice::span<const double> itemsToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemsToSend.size(), rankReceiver);
  MPI_Request request;
  MPI_Isend(const_cast<double *>(itemsToSend.data()),
            itemsToSend.size(), MPI_DOUBLE, rankReceiver, 0, _commState->comm, &request);
  return PtrRequest(new MPIRequest(request));
}

void MPIIntraComm::send(double itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);
  MPI_Send(&itemToSend, 1, MPI_DOUBLE, rankReceiver, 0, _commState->comm);
}

PtrRequest MPIIntraComm::aSend(const double &itemToSend, Rank rankReceiver)
{
  return aSend(precice::refToSpan<const double>(itemToSend), rankReceiver);
}

void MPIIntraComm::send(int itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);
  MPI_Send(&itemToSend, 1, MPI_INT, rankReceiver, 0, _commState->comm);
}

PtrRequest MPIIntraComm::aSend(const int &itemToSend, Rank rankReceiver)
{
  return aSend(precice::refToSpan<const int>(itemToSend), rankReceiver);
}

void MPIIntraComm::send(bool itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);
  MPI_Send(&itemToSend, 1, MPI_BOOL_INTRA, rankReceiver, 0, _commState->comm);
}

PtrRequest MPIIntraComm::aSend(const bool &itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE();
  MPI_Request request;
  MPI_Isend(const_cast<bool *>(&itemToSend),
            1, MPI_BOOL_INTRA, rankReceiver, 0, _commState->comm, &request);
  return PtrRequest(new MPIRequest(request));
}

void MPIIntraComm::receive(std::string &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  int        length;
  MPI_Status status;
  MPI_Probe(rankSender, 0, _commState->comm, &status);
  MPI_Get_count(&status, MPI_CHAR, &length);
  itemToReceive = std::string(length, '\0');
  MPI_Recv(const_cast<char *>(itemToReceive.data()),
           length, MPI_CHAR, rankSender, 0, _commState->comm, MPI_STATUS_IGNORE);
}

void MPIIntraComm::receive(precice::span<int> itemsToReceive, Rank rankSender)
{
  PRECICE_TRACE(itemsToReceive.size());
  MPI_Recv(itemsToReceive.data(),
           itemsToReceive.size(), MPI_INT, rankSender, 0, _commState->comm, MPI_STATUS_IGNORE);
}

void MPIIntraComm::receive(precice::span<double> itemsToReceive, Rank rankSender)
{
  PRECICE_TRACE(itemsToReceive.size());
  MPI_Recv(itemsToReceive.data(),
           itemsToReceive.size(), MPI_DOUBLE, rankSender, 0, _commState->comm, MPI_STATUS_IGNORE);
}

PtrRequest MPIIntraComm::aReceive(precice::span<double> itemsToReceive, int rankSender)
{
  PRECICE_TRACE(itemsToReceive.size());
  MPI_Request request;
  MPI_Irecv(itemsToReceive.data(),
            itemsToReceive.size(), MPI_DOUBLE, rankSender, 0, _commState->comm, &request);
  return PtrRequest(new MPIRequest(request));
}

void MPIIntraComm::receive(double &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  MPI_Recv(&itemToReceive, 1, MPI_DOUBLE, rankSender, 0, _commState->comm, MPI_STATUS_IGNORE);
}

PtrRequest MPIIntraComm::aReceive(double &itemToReceive, Rank rankSender)
{
  return aReceive(precice::refToSpan<double>(itemToReceive), rankSender);
}

void MPIIntraComm::receive(int &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  MPI_Recv(&itemToReceive, 1, MPI_INT, rankSender, 0, _commState->comm, MPI_STATUS_IGNORE);
}

PtrRequest MPIIntraComm::aReceive(int &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  MPI_Request request;
  MPI_Irecv(&itemToReceive, 1, MPI_INT, rankSender, 0, _commState->comm, &request);
  return PtrRequest(new MPIRequest(request));
}

void MPIIntraComm::receive(bool &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  MPI_Recv(&itemToReceive, 1, MPI_BOOL_INTRA, rankSender, 0, _commState->comm, MPI_STATUS_IGNORE);
}

PtrRequest MPIIntraComm::aReceive(bool &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  MPI_Request request;
  MPI_Irecv(&itemToReceive, 1, MPI_BOOL_INTRA, rankSender, 0, _commState->comm, &request);
  return PtrRequest(new MPIRequest(request));
}

} // namespace precice::com

#endif // not PRECICE_NO_MPI
