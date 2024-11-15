#pragma once
#include "time/Sample.hpp"

namespace precice::time {

/// @brief Stample containing timestampled Sample
struct Stample {
  double timestamp;
  Sample sample;
};

} // namespace precice::time
