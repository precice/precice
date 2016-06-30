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
