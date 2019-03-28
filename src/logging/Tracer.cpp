#include "Tracer.hpp"
#include <boost/log/attributes/mutable_constant.hpp>

namespace precice {
namespace logging {

Tracer::Tracer
(
  Logger &log,
  LogLocation loc
  )
  :
  _log(log),
  _loc(std::move(loc))
{}

Tracer::~Tracer()
{
  _log.trace(_loc, std::string{"Leaving "}.append(_loc.func));
}

}} // namespace precice,logging
