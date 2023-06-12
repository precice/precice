#include "Ginkgo.hpp"

#include "logging/Logger.hpp"

#ifndef PRECICE_NO_GINKGO
#include <Kokkos_Core.hpp>

namespace precice::utils {

bool Ginkgo::needs_finalize = false;

void Ginkgo::initialize(int *argc, char ***argv)
{
  if (!Kokkos::is_initialized()) {
    Kokkos::initialize(*argc, *argv);
    needs_finalize = true;
  }
}

void Ginkgo::finalize()
{
  if (needs_finalize) {
    Kokkos::finalize();
    needs_finalize = false;
  }
}

}

#endif