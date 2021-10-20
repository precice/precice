#include "Tracer.hpp"
#include <string>
#include <utility>
#include "logging/Logger.hpp"

namespace precice {
namespace logging {

Tracer::Tracer(
    Logger &    log,
    LogLocation loc)
    : _log(log),
      _loc(std::move(loc))
{
}

Tracer::~Tracer()
{
  _log.trace(_loc, std::string{"Leaving "}.append(_loc.func));
}

} // namespace logging
} // namespace precice
