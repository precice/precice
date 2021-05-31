#pragma once

#include <string>

namespace precice {
namespace cplscheme {
namespace constants {

const std::string &actionWriteIterationCheckpoint();

const std::string &actionReadIterationCheckpoint();

const std::string &actionWriteInitialData();

enum TimesteppingMethod {
  FIXED_TIME_WINDOW_SIZE,
  FIRST_PARTICIPANT_SETS_TIME_WINDOW_SIZE
};

} // namespace constants
} // namespace cplscheme
} // namespace precice
