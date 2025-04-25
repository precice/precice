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
    weInitialized = true;
  }
}

void Ginkgo::initialize(int nThreads, int deviceId)
{
  // We initialize Ginkgo internally through Kokkos
  if (!Kokkos::is_initialized() && !Kokkos::is_finalized()) {
    const int autoDeviceID = -1;
    // Strategy to select a device automatically from the GPUs available for execution. Must be either "mpi_rank" for round-robin assignment based on the local MPI rank or "random".
    // If kokkos was compiled with GPU as one backend, the GPU will always be initialized
    if (deviceId == autoDeviceID)
      Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(nThreads).set_map_device_id_by("mpi_rank").set_disable_warnings(true).set_print_configuration(false));
    else {
      Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(nThreads).set_device_id(deviceId).set_disable_warnings(true).set_print_configuration(false));
    }
    weInitialized = true;
  }
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
