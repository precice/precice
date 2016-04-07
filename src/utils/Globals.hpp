#pragma once

#include "tarch/Assertions.h"


#include "logging/LogMacros.hpp"
#include "LogMacros.hpp"

#include <string>

namespace precice {
namespace utils {

/**
 * @brief Holds global information.
 */
class Globals
{
public:

  static const std::string & getPathToSources ();

  static void setPathToSources ( const std::string & path );

private:

  static std::string _pathToSources;
};

}} // namespace precice, utils
