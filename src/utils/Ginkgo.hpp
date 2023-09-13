#pragma once

namespace precice {
namespace utils {

class Ginkgo {
public:
  static void initialize(int *argc, char ***argv);

  static void initialize(int nThreads, int deviceId);

  static void finalize();
};

} // namespace utils
} // namespace precice
