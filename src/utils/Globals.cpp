// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Globals.hpp"
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

