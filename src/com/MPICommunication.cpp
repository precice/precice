#ifndef PRECICE_NO_MPI

#include <cstddef>
#include <ostream>

#include <numeric>
#include "com/MPICommunication.hpp"
#include "com/MPIRequest.hpp"
#include "logging/LogMacros.hpp"
#include "precice/impl/Types.hpp"
#include "utils/span_tools.hpp"
#include "utils/IntraComm.hpp"

template <size_t>
struct MPI_Select_unsigned_integer_datatype;

template <>
struct MPI_Select_unsigned_integer_datatype<1> {
  static MPI_Datatype datatype;
};
MPI_Datatype MPI_Select_unsigned_integer_datatype<1>::datatype = MPI_UNSIGNED_CHAR;

template <>
struct MPI_Select_unsigned_integer_datatype<2> {
  static MPI_Datatype datatype;
};
MPI_Datatype MPI_Select_unsigned_integer_datatype<2>::datatype = MPI_UNSIGNED_SHORT;

template <>
struct MPI_Select_unsigned_integer_datatype<4> {
  static MPI_Datatype datatype;
};
MPI_Datatype MPI_Select_unsigned_integer_datatype<4>::datatype = MPI_UNSIGNED;

template <>
struct MPI_Select_unsigned_integer_datatype<8> {
  static MPI_Datatype datatype;
};
MPI_Datatype MPI_Select_unsigned_integer_datatype<8>::datatype = MPI_UNSIGNED_LONG;

#define MPI_BOOL MPI_Select_unsigned_integer_datatype<sizeof(bool)>::datatype

namespace precice::com {
MPICommunication::MPICommunication() = default;

void MPICommunication::send(std::string const &itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);
  rankReceiver = adjustRank(rankReceiver);
  PRECICE_DEBUG("Message: " + itemToSend);
  MPI_Send(const_cast<char *>(itemToSend.c_str()),
           itemToSend.size(),
           MPI_CHAR,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

void MPICommunication::send(precice::span<const int> itemsToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemsToSend.size());
  rankReceiver = adjustRank(rankReceiver);
  MPI_Send(const_cast<int *>(itemsToSend.data()),
           itemsToSend.size(),
           MPI_INT,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

PtrRequest MPICommunication::aSend(precice::span<const int> itemsToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemsToSend.size());
  rankReceiver = adjustRank(rankReceiver);

  MPI_Request request;
  MPI_Isend(const_cast<int *>(itemsToSend.data()),
            itemsToSend.size(),
            MPI_INT,
            rank(rankReceiver),
            0,
            communicator(rankReceiver),
            &request);

  return PtrRequest(new MPIRequest(request));
}

void MPICommunication::send(precice::span<const double> itemsToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemsToSend.size());
  rankReceiver = adjustRank(rankReceiver);
  MPI_Send(const_cast<double *>(itemsToSend.data()),
           itemsToSend.size(),
           MPI_DOUBLE,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

PtrRequest MPICommunication::aSend(precice::span<const double> itemsToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemsToSend.size(), rankReceiver);
  rankReceiver = adjustRank(rankReceiver);

  MPI_Request request;
  MPI_Isend(const_cast<double *>(itemsToSend.data()),
            itemsToSend.size(),
            MPI_DOUBLE,
            rank(rankReceiver),
            0,
            communicator(rankReceiver),
            &request);

  return PtrRequest(new MPIRequest(request));
}

void MPICommunication::send(double itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);
  rankReceiver = adjustRank(rankReceiver);
  MPI_Send(&itemToSend,
           1,
           MPI_DOUBLE,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

PtrRequest MPICommunication::aSend(const double &itemToSend, Rank rankReceiver)
{
  return aSend(precice::refToSpan<const double>(itemToSend), rankReceiver);
}

void MPICommunication::send(int itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);
  rankReceiver = adjustRank(rankReceiver);
  MPI_Send(&itemToSend,
           1,
           MPI_INT,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

PtrRequest MPICommunication::aSend(const int &itemToSend, Rank rankReceiver)
{
  return aSend(precice::refToSpan<const int>(itemToSend), rankReceiver);
}

void MPICommunication::send(bool itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);
  rankReceiver = adjustRank(rankReceiver);
  MPI_Send(&itemToSend,
           1,
           MPI_BOOL,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

PtrRequest MPICommunication::aSend(const bool &itemToSend, Rank rankReceiver)
{
  PRECICE_TRACE();
  rankReceiver = adjustRank(rankReceiver);

  MPI_Request request;
  MPI_Isend(const_cast<bool *>(&itemToSend),
            1,
            MPI_BOOL,
            rank(rankReceiver),
            0,
            communicator(rankReceiver),
            &request);

  return PtrRequest(new MPIRequest(request));
}

void MPICommunication::receive(std::string &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(itemToReceive, rankSender);
  rankSender = adjustRank(rankSender);
  int        length;
  MPI_Status status;
  MPI_Probe(rank(rankSender), 0, communicator(rankSender), &status);
  MPI_Get_count(&status, MPI_CHAR, &length);
  PRECICE_DEBUG("Stringlength = {}", length);
  itemToReceive = std::string(length, '\0');
  MPI_Recv(const_cast<char *>(itemToReceive.data()),
           length,
           MPI_CHAR,
           rank(rankSender),
           0,
           communicator(rankSender),
           MPI_STATUS_IGNORE);
  PRECICE_DEBUG("Received \"{}\" from rank {}", itemToReceive, rankSender);
}

void MPICommunication::receive(precice::span<int> itemsToReceive, Rank rankSender)
{
  PRECICE_TRACE(itemsToReceive.size());
  rankSender = adjustRank(rankSender);

  MPI_Status status;
  MPI_Recv(itemsToReceive.data(),
           itemsToReceive.size(),
           MPI_INT,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
}

void MPICommunication::receive(precice::span<double> itemsToReceive, Rank rankSender)
{
  PRECICE_TRACE(itemsToReceive.size());
  rankSender = adjustRank(rankSender);

  MPI_Status status;
  MPI_Recv(itemsToReceive.data(),
           itemsToReceive.size(),
           MPI_DOUBLE,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
}

PtrRequest MPICommunication::aReceive(precice::span<double> itemsToReceive, Rank rankSender)
{
  PRECICE_TRACE(itemsToReceive.size());
  rankSender = adjustRank(rankSender);

  MPI_Request request;
  MPI_Irecv(itemsToReceive.data(),
            itemsToReceive.size(),
            MPI_DOUBLE,
            rank(rankSender),
            0,
            communicator(rankSender),
            &request);

  return PtrRequest(new MPIRequest(request));
}

void MPICommunication::receive(double &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  rankSender = adjustRank(rankSender);

  MPI_Status status;
  MPI_Recv(&itemToReceive,
           1,
           MPI_DOUBLE,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
  PRECICE_DEBUG("Received {} from rank {}", itemToReceive, rankSender);
}

PtrRequest MPICommunication::aReceive(double &itemToReceive, Rank rankSender)
{
  return aReceive(precice::refToSpan<double>(itemToReceive), rankSender);
}

void MPICommunication::receive(int &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  rankSender = adjustRank(rankSender);

  MPI_Status status;
  MPI_Recv(&itemToReceive,
           1,
           MPI_INT,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
  PRECICE_DEBUG("Received {} from rank {}", itemToReceive, rankSender);
}

PtrRequest MPICommunication::aReceive(int &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  rankSender = adjustRank(rankSender);

  MPI_Request request;
  MPI_Irecv(&itemToReceive,
            1,
            MPI_INT,
            rank(rankSender),
            0,
            communicator(rankSender),
            &request);

  return PtrRequest(new MPIRequest(request));
}

void MPICommunication::receive(bool &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  rankSender = adjustRank(rankSender);

  MPI_Status status;
  MPI_Recv(&itemToReceive,
           1,
           MPI_BOOL,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
  PRECICE_DEBUG("Received {} from rank {}", itemToReceive, rankSender);
}

PtrRequest MPICommunication::aReceive(bool &itemToReceive, Rank rankSender)
{
  PRECICE_TRACE(rankSender);
  rankSender = adjustRank(rankSender);

  MPI_Request request;
  MPI_Irecv(&itemToReceive,
            1,
            MPI_BOOL,
            rank(rankSender),
            0,
            communicator(rankSender),
            &request);

  return PtrRequest(new MPIRequest(request));
}

void MPICommunication::gather(int itemToSend, std::vector<int> &itemsToReceive)
{
  PRECICE_TRACE(itemToSend);

  Rank rootRank = adjustRank(0);

  MPI_Gather(&itemToSend,
             1,
             MPI_INT,
             itemsToReceive.data(),
             1,
             MPI_INT,
             rootRank,
             communicator(rootRank) // TODO: We cannot just assume that the communicator is the same for all!
  );
}

void MPICommunication::gather(span<const int> itemToSend, std::vector<std::vector<int>>& itemsToReceive, const std::vector<int>& recvcounts)
{
  PRECICE_TRACE(itemToSend.size(), recvcounts);

  Rank rootRank = adjustRank(0);
  bool isPrimary = utils::IntraComm::isPrimary();

  std::vector<int> flatBuffer;
  std::vector<int> displs;

  if (isPrimary) {
    displs.resize(recvcounts.size());
    int totalSize = std::accumulate(recvcounts.begin(), recvcounts.end(), 0);

    std::exclusive_scan(
        recvcounts.begin(),
        recvcounts.end(),
        displs.begin(),
        0);
    flatBuffer.resize(totalSize);
  }

  PRECICE_WARN("This gather only works if all ranks are in the same communicator. It will likely crash when using mpi-multiple-ports (MPIPortsCommunication instead of MPISinglePortsCommunication)");

  MPI_Gatherv(itemToSend.data(),
              itemToSend.size(),
              MPI_INT,
              isPrimary ? flatBuffer.data() : nullptr,
              isPrimary ? recvcounts.data() : nullptr,
              isPrimary ? displs.data() : nullptr,
              MPI_INT,
              rootRank,
              communicator(rootRank) // TODO: We cannot just assume that the communicator is the same for all!
  );

  if (isPrimary) {
    for (int i = 0; i < recvcounts.size(); i++) {
      itemsToReceive[i] = std::vector<int>(
        flatBuffer.begin() + displs[i],
        flatBuffer.begin() + displs[i] + recvcounts[i]
      );
    }
  }
}

} // namespace precice::com

#endif // not PRECICE_NO_MPI
