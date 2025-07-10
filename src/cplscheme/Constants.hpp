#pragma once

#include <string>

namespace precice::cplscheme::constants {

enum TimesteppingMethod {
  FIXED_TIME_WINDOW_SIZE,
  FIRST_PARTICIPANT_SETS_TIME_WINDOW_SIZE
};

} // namespace precice::cplscheme::constants
