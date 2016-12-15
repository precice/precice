#pragma once

#include "MasterSlave.hpp"

#include "utils/assertion.hpp"
#include "utils/prettyprint.hpp"

#include "logging/LogMacros.hpp"

#include <string>

namespace precice {
namespace utils {

/// Holds global information.
class Globals
{
public:

  static const std::string & getPathToSources ();

  static void setPathToSources ( const std::string & path );

private:

  static std::string _pathToSources;
};

}} // namespace precice, utils
