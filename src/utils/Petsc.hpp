#pragma once

#include "logging/Logger.hpp"

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
  
  /// Whether we have initialized Petsc or if it was initialized by an application calling us.
  static bool weInitialized;
  
  static logging::Logger _log;
};


}} // namespace precice, utils

