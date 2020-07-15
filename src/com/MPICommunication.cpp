#ifndef PRECICE_NO_MPI

#include "MPICommunication.hpp"
#include <ostream>
#include <stddef.h>
#include "MPIRequest.hpp"
#include "logging/LogMacros.hpp"

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

namespace precice {
namespace com {
MPICommunication::MPICommunication()
{
}

void MPICommunication::send(std::string const &itemToSend, int rankReceiver)
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

void MPICommunication::send(const int *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size);
  rankReceiver = adjustRank(rankReceiver);
  MPI_Send(const_cast<int *>(itemsToSend),
           size,
           MPI_INT,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

PtrRequest MPICommunication::aSend(const int *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size);
  rankReceiver = adjustRank(rankReceiver);

  MPI_Request request;
  MPI_Isend(const_cast<int *>(itemsToSend),
            size,
            MPI_INT,
            rank(rankReceiver),
            0,
            communicator(rankReceiver),
            &request);

  return PtrRequest(new MPIRequest(request));
}

void MPICommunication::send(const double *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size);
  rankReceiver = adjustRank(rankReceiver);
  MPI_Send(const_cast<double *>(itemsToSend),
           size,
           MPI_DOUBLE,
           rank(rankReceiver),
           0,
           communicator(rankReceiver));
}

PtrRequest MPICommunication::aSend(const double *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size, rankReceiver);
  rankReceiver = adjustRank(rankReceiver);

  MPI_Request request;
  MPI_Isend(const_cast<double *>(itemsToSend),
            size,
            MPI_DOUBLE,
            rank(rankReceiver),
            0,
            communicator(rankReceiver),
            &request);

  return PtrRequest(new MPIRequest(request));
}

PtrRequest MPICommunication::aSend(std::vector<double> const &itemsToSend, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver, itemsToSend.size(), itemsToSend);
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

void MPICommunication::send(double itemToSend, int rankReceiver)
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

PtrRequest MPICommunication::aSend(const double &itemToSend, int rankReceiver)
{
  return aSend(&itemToSend, 1, rankReceiver);
}

void MPICommunication::send(int itemToSend, int rankReceiver)
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

PtrRequest MPICommunication::aSend(const int &itemToSend, int rankReceiver)
{
  return aSend(&itemToSend, 1, rankReceiver);
}

void MPICommunication::send(bool itemToSend, int rankReceiver)
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

PtrRequest MPICommunication::aSend(const bool &itemToSend, int rankReceiver)
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

void MPICommunication::receive(std::string &itemToReceive, int rankSender)
{
  PRECICE_TRACE(itemToReceive, rankSender);
  rankSender = adjustRank(rankSender);
  int        length;
  MPI_Status status;
  MPI_Probe(rank(rankSender), 0, communicator(rankSender), &status);
  MPI_Get_count(&status, MPI_CHAR, &length);
  PRECICE_DEBUG("Stringlength = " << length);
  itemToReceive = std::string(length, '\0');
  MPI_Recv(const_cast<char *>(itemToReceive.data()),
           length,
           MPI_CHAR,
           rank(rankSender),
           0,
           communicator(rankSender),
           MPI_STATUS_IGNORE);
  PRECICE_DEBUG("Received \"" << itemToReceive << "\" from rank " << rankSender);
}

void MPICommunication::receive(int *itemsToReceive, int size, int rankSender)
{
  PRECICE_TRACE(size);
  rankSender = adjustRank(rankSender);

  MPI_Status status;
  MPI_Recv(itemsToReceive,
           size,
           MPI_INT,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
}

void MPICommunication::receive(double *itemsToReceive, int size, int rankSender)
{
  PRECICE_TRACE(size);
  rankSender = adjustRank(rankSender);

  MPI_Status status;
  MPI_Recv(itemsToReceive,
           size,
           MPI_DOUBLE,
           rank(rankSender),
           0,
           communicator(rankSender),
           &status);
}

PtrRequest MPICommunication::aReceive(double *itemsToReceive, int size, int rankSender)
{
  PRECICE_TRACE(size);
  rankSender = adjustRank(rankSender);

  MPI_Request request;
  MPI_Irecv(itemsToReceive,
            size,
            MPI_DOUBLE,
            rank(rankSender),
            0,
            communicator(rankSender),
            &request);

  return PtrRequest(new MPIRequest(request));
}

PtrRequest MPICommunication::aReceive(std::vector<double> &itemsToReceive, int rankSender)
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

void MPICommunication::receive(double &itemToReceive, int rankSender)
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
  PRECICE_DEBUG("Received " << itemToReceive << " from rank " << rankSender);
}

PtrRequest MPICommunication::aReceive(double &itemToReceive, int rankSender)
{
  return aReceive(&itemToReceive, 1, rankSender);
}

void MPICommunication::receive(int &itemToReceive, int rankSender)
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
  PRECICE_DEBUG("Received " << itemToReceive << " from rank " << rankSender);
}

PtrRequest MPICommunication::aReceive(int &itemToReceive, int rankSender)
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

void MPICommunication::receive(bool &itemToReceive, int rankSender)
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
  PRECICE_DEBUG("Received " << itemToReceive << " from rank " << rankSender);
}

PtrRequest MPICommunication::aReceive(bool &itemToReceive, int rankSender)
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

void MPICommunication::send(std::vector<int> const &v, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver);
  rankReceiver = adjustRank(rankReceiver);
  MPI_Send(const_cast<int *>(v.data()), v.size(), MPI_INT,
           rank(rankReceiver), 0, communicator(rankReceiver));
}

void MPICommunication::receive(std::vector<int> &v, int rankSender)
{
  PRECICE_TRACE(rankSender);
  rankSender        = adjustRank(rankSender);
  int        length = -1;
  MPI_Status status;
  MPI_Probe(rank(rankSender), 0, communicator(rankSender), &status);
  MPI_Get_count(&status, MPI_INT, &length);
  v.resize(length);
  MPI_Recv(v.data(), length, MPI_INT, rank(rankSender),
           0, communicator(rankSender), MPI_STATUS_IGNORE);
}

void MPICommunication::send(std::vector<double> const &v, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver);
  rankReceiver = adjustRank(rankReceiver);
  MPI_Send(const_cast<double *>(v.data()), v.size(), MPI_DOUBLE,
           rank(rankReceiver), 0, communicator(rankReceiver));
}

void MPICommunication::receive(std::vector<double> &v, int rankSender)
{
  PRECICE_TRACE(rankSender);
  rankSender        = adjustRank(rankSender);
  int        length = -1;
  MPI_Status status;
  MPI_Probe(rank(rankSender), 0, communicator(rankSender), &status);
  MPI_Get_count(&status, MPI_DOUBLE, &length);
  v.resize(length, '\0');
  MPI_Recv(v.data(), length, MPI_DOUBLE, rank(rankSender),
           0, communicator(rankSender), MPI_STATUS_IGNORE);
}

} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
