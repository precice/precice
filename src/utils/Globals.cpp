#include "utils/Globals.hpp"
#include <cstdlib>
#include <string>

namespace precice {
namespace utils {

std::string getPathToSources()
{
  logging::Logger _log("utils");
  char* preciceRoot = std::getenv("PRECICE_ROOT");
  CHECK(preciceRoot != nullptr,
        "Environment variable PRECICE_ROOT has not been set. Please set it to the precice directory.");
  return std::string(preciceRoot) + "/src";
}

}} // namespace precice, utils

