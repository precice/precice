#include <string>
#include <utility>

#include "logging/Tracer.hpp"

namespace precice::logging {

Tracer::Tracer(
    Logger     &log,
    LogLocation loc)
    : _log(log),
      _loc(std::move(loc))
{
}

Tracer::~Tracer()
{
  _log.trace(_loc, std::string{"Leaving "}.append(_loc.func));
}

} // namespace precice::logging
