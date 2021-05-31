#pragma once

/**
 * MPI Mock
 *
 * Inlines because of getting multiple definitons from the compiler otherwise.
 */

using MPI_Comm     = std::nullptr_t;
using MPI_Op       = std::nullptr_t;
using MPI_Datatype = std::nullptr_t;

static MPI_Comm MPI_COMM_WORLD = nullptr;

const MPI_Datatype MPI_LONG = nullptr;
const MPI_Op       MPI_MIN  = nullptr;
const MPI_Op       MPI_MAX  = nullptr;

inline int MPI_Barrier(MPI_Comm comm)
{
  return 0;
}

inline int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
  *rank = 0;
  return 0;
}

inline int MPI_Comm_size(MPI_Comm comm, int *size)
{
  *size = 0;
  return 0;
}

template <class T>
inline int MPI_Allreduce(const T *sendbuf, T *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  std::copy(sendbuf, sendbuf + count, recvbuf);
  return 0;
}

template <class T>
inline int MPI_Reduce(const T *sendbuf, T *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
  std::copy(sendbuf, sendbuf + count, recvbuf);
  return 0;
}
