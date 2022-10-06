#pragma once
#ifndef PRECICE_NO_MPI

#include <mpi.h>
#include "com/Request.hpp"

namespace precice {
namespace com {
class MPIRequest : public Request {
public:
  explicit MPIRequest(MPI_Request request);

  bool test() override;

  void wait() override;

  boost::system::error_code errorCode() const noexcept override
  {
    return _ec;
  }

private:
  boost::system::error_code _ec;
  MPI_Request               _request;
};
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
