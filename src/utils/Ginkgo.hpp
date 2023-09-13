#pragma once

namespace precice {
namespace utils {

class Ginkgo {
public:
  static void initialize(int *argc, char ***argv);

  static void finalize();
};

} // namespace utils
} // namespace precice
