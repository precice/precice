#include "Testing.hpp"
#include "logging/LogMacros.hpp"
#include <cstdlib>

namespace precice {
namespace testing {

std::string getPathToSources()
{
  precice::logging::Logger _log("testing");
  char* preciceRoot = std::getenv("PRECICE_ROOT");
  P_CHECK(preciceRoot != nullptr,
        "Environment variable PRECICE_ROOT has not been set. Please set it to the precice directory.");
  return std::string(preciceRoot) + "/src";
}


}
}
