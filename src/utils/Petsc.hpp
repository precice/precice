#pragma once

#include "tarch/logging/Log.h"

namespace precice {
namespace utils {

/// Utility class for managing PETSc operations.
class Petsc
{
public:



  /**
   * @brief Initializes the Petsc environment.
   *
   * @param[in] argc Parameter count
   * @param[in] argc Parameter values, is passed to PetscInitialize
   */
  static void initialize (
    int*               argc,
    char***            argv);

  /// Finalizes Petsc environment.
  static void finalize();


private:

  static tarch::logging::Log _log;
};


}} // namespace precice, utils

