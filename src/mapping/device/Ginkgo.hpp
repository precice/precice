#pragma once

namespace precice {
namespace device {

class Ginkgo {
public:
  static void initialize(int *argc, char ***argv);

  static void initialize(int nThreads, int deviceId);

  static void finalize();

private:
  /// Whether we have initialized Ginkgo or if it was initialized by an application calling us.
  static bool weInitialized;
};

} // namespace device
} // namespace precice
