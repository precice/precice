#pragma once

#include "precice/exceptions.hpp"
#include "utils/fmt.hpp"

#define PRECICE_LOG_LOCATION     \
  precice::logging::LogLocation  \
  {                              \
    __FILE__, __LINE__, __func__ \
  }

#define PRECICE_WARN(...) _log.warning(PRECICE_LOG_LOCATION, precice::utils::format_or_error(__VA_ARGS__))

#define PRECICE_INFO(...) _log.info(PRECICE_LOG_LOCATION, precice::utils::format_or_error(__VA_ARGS__))

#define PRECICE_ERROR(...)                                                   \
  do {                                                                       \
    auto preciceErrorMessage = precice::utils::format_or_error(__VA_ARGS__); \
    _log.error(PRECICE_LOG_LOCATION, preciceErrorMessage);                   \
    throw precice::Error{preciceErrorMessage};                               \
  } while (false)

#define PRECICE_CHECK(check, ...) \
  if (!(check)) {                 \
    PRECICE_ERROR(__VA_ARGS__);   \
  }

// Debug logging is disabled in release (NDEBUG) builds by default.
// To enable it anyhow, enable the CMake option PRECICE_RELEASE_WITH_DEBUG_LOG.

#if defined(NDEBUG) && !defined(PRECICE_RELEASE_WITH_DEBUG_LOG)
#define PRECICE_NO_DEBUG_LOG
#endif

#ifdef PRECICE_NO_DEBUG_LOG

#define PRECICE_DEBUG(...) \
  {                        \
  }
#define PRECICE_TRACE(...) \
  {                        \
  }

#else // PRECICE_NO_DEBUG_LOG

#include "utils/ArgumentFormatter.hpp"

#define PRECICE_DEBUG(...) _log.debug(PRECICE_LOG_LOCATION, precice::utils::format_or_error(__VA_ARGS__))

#endif // ! PRECICE_NO_DEBUG_LOG

// Trace logging is disabled in release (NDEBUG) builds by default.
// To enable it anyhow, enable the CMake option PRECICE_RELEASE_WITH_TRACE_LOG.

#if defined(NDEBUG) && !defined(PRECICE_RELEASE_WITH_TRACE_LOG)
#define PRECICE_NO_TRACE_LOG
#endif

#ifdef PRECICE_NO_TRACE_LOG

#define PRECICE_TRACE(...) \
  {                        \
  }

#else // PRECICE_NO_TRACE_LOG

#include "logging/Tracer.hpp"

// Do not put do {...} while (false) here, it will destroy the _tracer_ right after creation
#define PRECICE_TRACE(...)                                       \
  precice::logging::Tracer _tracer_(_log, PRECICE_LOG_LOCATION); \
  _log.trace(PRECICE_LOG_LOCATION, std::string{"Entering "} + __func__ + PRECICE_LOG_ARGUMENTS(__VA_ARGS__))

#endif // ! PRECICE_NO_TRACE_LOG
