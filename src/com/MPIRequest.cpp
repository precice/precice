#ifndef PRECICE_NO_MPI

#include "com/MPIRequest.hpp"

namespace precice::com {
MPIRequest::MPIRequest(MPI_Request request)
    : _request(request)
{
}

bool MPIRequest::test()
{
  int        complete = 0;
  MPI_Status status;
  MPI_Test(&_request, &complete, &status);

  if (complete) {
    int cancelled = 0;
    MPI_Test_cancelled(&status, &cancelled);

    _ec = boost::system::errc::make_error_code(
        cancelled ? boost::system::errc::operation_canceled : boost::system::errc::success);
  }

  return complete;
}

void MPIRequest::wait()
{
  MPI_Status status;
  MPI_Wait(&_request, &status);

  int cancelled = 0;
  MPI_Test_cancelled(&status, &cancelled);

  _ec = boost::system::errc::make_error_code(
      cancelled ? boost::system::errc::operation_canceled : boost::system::errc::success);
}
} // namespace precice::com

#endif // not PRECICE_NO_MPI
