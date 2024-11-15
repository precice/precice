#pragma once

#include <string>

namespace precice {
namespace cplscheme {
namespace constants {

enum TimesteppingMethod {
  FIXED_TIME_WINDOW_SIZE,
  FIRST_PARTICIPANT_SETS_TIME_WINDOW_SIZE
};

} // namespace constants
} // namespace cplscheme
} // namespace precice
