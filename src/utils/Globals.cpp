#include "utils/Globals.hpp"

namespace precice {
namespace utils {

std::string Globals:: _pathToSources = "";

const std::string & Globals:: getPathToSources ()
{
  assertion ( not _pathToSources.empty() );
  return _pathToSources;
}

void Globals:: setPathToSources
(
  const std::string & path )
{
  assertion ( not path.empty() );
  _pathToSources = path;
}

}} // namespace precice, utils

