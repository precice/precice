#include "Ginkgo.hpp"

#include "logging/Logger.hpp"

#ifndef PRECICE_NO_GINKGO
#include <Kokkos_Core.hpp>

namespace precice::utils {

// bool Ginkgo::needs_finalize = false;

void Ginkgo::initialize(int *argc, char ***argv)
{
  if (!Kokkos::is_initialized() && !Kokkos::is_finalized()) {
    Kokkos::initialize(*argc, *argv);
  }
}

void Ginkgo::finalize()
{
  if (Kokkos::is_initialized()) {
    Kokkos::finalize();
  }
}

} // namespace precice::utils

#endif
