#pragma once

#include <memory>
#include <string_view>

namespace precice::logging {

/// Struct used to capture the original location of a log request
struct LogLocation {
  const char *file;
  int         line;
  const char *func;
};

/// This class provides a lightweight logger.
class Logger {
public:
  /** Creates a logger for a given module.
   * @param[in] the name of the module
   */
  explicit Logger(std::string_view module);

  Logger(const Logger &other);
  Logger(Logger &&other) noexcept;
  Logger &operator=(Logger other);
  ~Logger();

  void swap(Logger &other) noexcept;

  ///@name Logging operations
  ///@{
  void error(LogLocation loc, std::string_view mess) noexcept;
  void warning(LogLocation loc, std::string_view mess) noexcept;
  void info(LogLocation loc, std::string_view mess) noexcept;
  void debug(LogLocation loc, std::string_view mess) noexcept;
  void trace(LogLocation loc, std::string_view mess) noexcept;
  ///@}

private:
  /// Forward declaration of the implementation of the logger
  class LoggerImpl;

  /// Pimpl to the logger implementation
  std::unique_ptr<LoggerImpl> _impl;
};

} // namespace precice::logging

// Include LogMacros here, because using it works only together with a Logger
#include "LogMacros.hpp"
