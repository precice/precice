#pragma once

#include "logging/Logger.hpp"

namespace precice::logging {

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

} // namespace precice::logging
