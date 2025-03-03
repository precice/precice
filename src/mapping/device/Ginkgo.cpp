#include "Ginkgo.hpp"

#include "logging/Logger.hpp"

#ifndef PRECICE_NO_GINKGO
#include <Kokkos_Core.hpp>

namespace precice::device {

bool Ginkgo::weInitialized = false;

void Ginkgo::initialize(int *argc, char ***argv)
{
  // We initialize Ginkgo internally through Kokkos
  if (!Kokkos::is_initialized() && !Kokkos::is_finalized()) {
    Kokkos::initialize(*argc, *argv);
  }
  weInitialized = true;
}

void Ginkgo::initialize(int nThreads, int deviceId)
{
  // We initialize Ginkgo internally through Kokkos
  if (!Kokkos::is_initialized() && !Kokkos::is_finalized()) {
    const int defaultDeviceID = -1;
    // don't use and initialize the GPU if we don't use it
    if (deviceId == defaultDeviceID)
      Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(nThreads).set_disable_warnings(true));
    else {
      Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(nThreads).set_device_id(deviceId).set_disable_warnings(true));
    }
  }
  weInitialized = true;
}

void Ginkgo::finalize()
{
  // we finalize internally through Kokkos as well
  if (weInitialized && Kokkos::is_initialized()) {
    Kokkos::finalize();
  }
}

} // namespace precice::device

#endif
