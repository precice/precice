#pragma once

#include <memory>
#include <string>

namespace precice {
namespace logging {

/// Struct used to capture the original location of a log request
struct LogLocation {
  const char *file;
  int         line;
  const char *func;
};

/// This class provides a leightweight logger.
class Logger {
public:
  /** Creates a logger for a given module.
   * @param[in] the name of the module 
   */
  explicit Logger(std::string module);

  Logger(const Logger &other);
  Logger(Logger &&other);
  Logger &operator=(Logger other);
  ~Logger();

  void swap(Logger &other) noexcept;

  ///@name Logging operations
  ///@{
  void error(LogLocation loc, const std::string &mess) noexcept;
  void warning(LogLocation loc, const std::string &mess) noexcept;
  void info(LogLocation loc, const std::string &mess) noexcept;
  void debug(LogLocation loc, const std::string &mess) noexcept;
  void trace(LogLocation loc, const std::string &mess) noexcept;
  ///@}

private:
  /// Forward declaration of the implementation of the logger
  class LoggerImpl;

  /// Pimpl to the logger implementation
  std::unique_ptr<LoggerImpl> _impl;
};

} // namespace logging
} // namespace precice

// Include LogMacros here, because using it works only together with a Logger
#include "LogMacros.hpp"
