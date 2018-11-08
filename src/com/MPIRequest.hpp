#pragma once
#ifndef PRECICE_NO_MPI

#include "Request.hpp"
#include <mpi.h>

namespace precice
{
namespace com
{
class MPIRequest : public Request
{
public:
  explicit MPIRequest(MPI_Request request);

  bool test() override;

  void wait() override;

private:
  MPI_Request _request;
};
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
