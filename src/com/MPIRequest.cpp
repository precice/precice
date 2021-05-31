#ifndef PRECICE_NO_MPI

#include "MPIRequest.hpp"

namespace precice {
namespace com {
MPIRequest::MPIRequest(MPI_Request request)
    : _request(request)
{
}

bool MPIRequest::test()
{
  int complete = 0;

  MPI_Test(&_request, &complete, MPI_STATUS_IGNORE);

  return complete;
}

void MPIRequest::wait()
{
  MPI_Wait(&_request, MPI_STATUS_IGNORE);
}
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
