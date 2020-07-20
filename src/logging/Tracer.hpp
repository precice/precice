#pragma once
#include "Logger.hpp"

namespace precice {
namespace logging {
class Logger;
struct LogLocation;

class Tracer {
public:
  Tracer(Logger &log, LogLocation loc);
  ~Tracer();

private:
  Logger &_log;

  LogLocation _loc;
};

} // namespace logging
} // namespace precice
