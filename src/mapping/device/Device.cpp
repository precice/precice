#include "Device.hpp"

#include "logging/Logger.hpp"

#if !defined(PRECICE_NO_GINKGO) || !defined(PRECICE_NO_KOKKOS_KERNELS)
#include <Kokkos_Core.hpp>

namespace precice::device {

bool Device::weInitialized = false;

void Device::initialize(int *argc, char ***argv)
{
  if (!Kokkos::is_initialized() && !Kokkos::is_finalized()) {
    Kokkos::initialize(*argc, *argv);
    weInitialized = true;
  }
}

void Device::initialize(int nThreads, int deviceId)
{
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

void Device::finalize()
{
  if (weInitialized && Kokkos::is_initialized()) {
    Kokkos::finalize();
  }
}

} // namespace precice::device

#endif
