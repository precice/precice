#pragma once

namespace precice {
namespace device {

class Ginkgo {
public:
  static void initialize(int *argc, char ***argv);

  static void initialize(int nThreads, int deviceId);

  static void finalize();
};

} // namespace device
} // namespace precice
