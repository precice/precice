#include "utils/Globals.hpp"
#include <cstdlib>
#include <string>

namespace precice {
namespace utils {

std::string getPathToSources()
{
  char* preciceRoot = std::getenv("PRECICE_ROOT");
  assertion(preciceRoot != nullptr,
            "Environment variable PRECICE_ROOT has not been set. Please set it to the precice directory.");
  return std::string(preciceRoot) + "/src";
}

}} // namespace precice, utils

