#include <ostream>
#ifndef PRECICE_NO_MPI

#pragma once

#include <mpi.h>
#include <ostream>

#include "utils/String.hpp"

namespace precice::utils {

/** Utility class to simplify use of MPI return codes
 *
 * Provides implicit conversion to bool and message printing.
 */
struct MPIResult {
  int code;

  MPIResult() = default;

  MPIResult(const MPIResult &other) = default;

  MPIResult(int res)
      : code(res) {}

  MPIResult &operator=(MPIResult other)
  {
    code = other.code;
    return *this;
  }

  operator bool() const
  {
    return code == MPI_SUCCESS;
  }

  std::string message() const
  {
    utils::StringMaker<MPI_MAX_ERROR_STRING> sm;

    int i;
    MPI_Error_string(code, sm.data(), &i);
    return sm.str();
  }
};

} // namespace precice::utils

#endif
