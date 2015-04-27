// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#include "MPIRequest.hpp"

namespace precice {
namespace com {
MPIRequest::MPIRequest(MPI_Request request) : _request(request) {
}

bool
MPIRequest::test() {
  int complete = 0;

  MPI_Test(&_request, &complete, MPI_STATUS_IGNORE);

  return complete;
}

void
MPIRequest::wait() {
  MPI_Wait(&_request, MPI_STATUS_IGNORE);
}
}
} // namespace precice, com

#endif // not PRECICE_NO_MPI
