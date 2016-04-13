#pragma once
#include "Logger.hpp"

namespace precice {
namespace logging {

class Tracer
{
public:  
  
  Tracer (Logger &log, std::string function, std::string file, long line);
  ~Tracer();

private:

  Logger _log;

  std::string _function;

  std::string _file;

  long _line;

};

}} // namespace precice, logging
