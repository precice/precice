#pragma once

#include <boost/stacktrace/stacktrace.hpp>
#include <string>

/// Returns a demangled stack backtrace of the caller function
inline auto getStacktrace() -> std::string
{
  std::ostringstream strm;
  try {
    strm << boost::stacktrace::stacktrace();
  } catch (const std::exception &e) {
    strm << "Stacktrace failed: " << e.what();
  }
  return strm.str();
}
