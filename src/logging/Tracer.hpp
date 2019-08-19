#pragma once
#include "Logger.hpp"

namespace precice {
namespace logging {

class Tracer
{
public:  
  
  Tracer (Logger &log, LogLocation loc);
  ~Tracer();

private:

  Logger _log;

  LogLocation _loc;
};

}} // namespace precice, logging
