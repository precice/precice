#pragma once

#include <string>

namespace precice {
namespace cplscheme {
namespace constants {

const std::string &actionWriteIterationCheckpoint();

const std::string &actionReadIterationCheckpoint();

const std::string &actionWriteInitialData();

enum TimesteppingMethod {
  FIXED_DT,
  FIRST_PARTICIPANT_SETS_DT
};

} // namespace constants
} // namespace cplscheme
} // namespace precice
