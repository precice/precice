#include "utils/Globals.hpp"
#include <cstdlib>

namespace precice {
namespace utils {

std::string getPathToSources()
{
  std::string root(std::getenv("PRECICE_ROOT"));
  assertion(not root.empty());
  return root + "/src";
}

}} // namespace precice, utils

