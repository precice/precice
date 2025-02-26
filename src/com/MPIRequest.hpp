#pragma once
#ifndef PRECICE_NO_MPI

#include <mpi.h>
#include "com/Request.hpp"

namespace precice::com {
class MPIRequest : public Request {
public:
  explicit MPIRequest(MPI_Request request);

  bool test() override;

  void wait() override;

private:
  MPI_Request _request;
};
} // namespace precice::com

#endif // not PRECICE_NO_MPI
