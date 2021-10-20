#include "testing/Testing.hpp"
#include <cstdlib>
#include <string>
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace testing {

std::string getPathToSources()
{
  precice::logging::Logger _log("testing");
  char *                   preciceRoot = std::getenv("PRECICE_ROOT");
  PRECICE_CHECK(preciceRoot != nullptr,
                "Environment variable PRECICE_ROOT is required to run the tests, but has not been set. Please set it to the root directory of the precice repository.");
  return std::string(preciceRoot) + "/src";
}

} // namespace testing
} // namespace precice
