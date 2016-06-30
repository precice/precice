#ifndef PRECICE_CPLSCHEME_CONSTANTS
#define PRECICE_CPLSCHEME_CONSTANTS

#include <string>

namespace precice {
namespace cplscheme {
namespace constants {

const std::string& actionWriteIterationCheckpoint();

const std::string& actionReadIterationCheckpoint();

const std::string& actionWriteInitialData();

enum TimesteppingMethod
{
  FIXED_DT,
  FIRST_PARTICIPANT_SETS_DT
};


}}} // namespace precice, cplscheme, constants

#endif // PRECICE_CPLSCHEME_CONSTANTS
