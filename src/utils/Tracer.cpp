// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Tracer.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace utils {

Tracer:: Tracer
(
  const tarch::logging::Log & log,
  const std::string & methodname,
  const std::string & stateString )
:
  _log ( log ),
  _methodname ( methodname ),
  _stateString ( stateString )
{
  std::string preciceMethodName (_methodname);
  preciceDebug ( "Entering " + stateString );
}

Tracer:: ~Tracer ()
{
  std::string preciceMethodName (_methodname);
  preciceDebug ( "Leaving" );
}

}} // namespace precice, utils
