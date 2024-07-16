#include "Ginkgo.hpp"

#include "logging/Logger.hpp"

#ifndef PRECICE_NO_GINKGO
#include <Kokkos_Core.hpp>

namespace precice::device {

// bool Ginkgo::needs_finalize = false;

void Ginkgo::initialize(int *argc, char ***argv)
{
  if (!Kokkos::is_initialized() && !Kokkos::is_finalized()) {
    Kokkos::initialize(*argc, *argv);
  }
}

void Ginkgo::initialize(int nThreads, int deviceId)
{
  if (!Kokkos::is_initialized() && !Kokkos::is_finalized()) {
    Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(nThreads).set_device_id(deviceId).set_disable_warnings(true));
  }
}

void Ginkgo::finalize()
{
  if (Kokkos::is_initialized()) {
    Kokkos::finalize();
  }
}

} // namespace precice::device

#endif
