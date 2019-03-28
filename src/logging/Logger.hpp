#pragma once

#include <string>
#include <memory>

namespace precice {
namespace logging {

struct LogLocation {
    const char * file;
    int line;
    const char * func;
};

class Logger {
public:
  explicit Logger(std::string module);
  Logger(const Logger& other);
  Logger& operator=(Logger other);
  ~Logger();

  void swap(Logger& other) noexcept;

  void error(LogLocation loc, const std::string& mess);
  void warning(LogLocation loc, const std::string& mess);
  void info(LogLocation loc, const std::string& mess);
  void debug(LogLocation loc, const std::string& mess);
  void trace(LogLocation loc, const std::string& mess);

private:
  class LoggerImpl;
  std::unique_ptr<LoggerImpl> _impl;
};

}} // namespace precice, logging

// Include LogMacros here, because using it works only together with a Logger
#include "LogMacros.hpp" 

