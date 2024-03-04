#pragma once

#include "utils/fmt.hpp"

#define PRECICE_LOG_LOCATION     \
  precice::logging::LogLocation  \
  {                              \
    __FILE__, __LINE__, __func__ \
  }

#define PRECICE_WARN(...) _log.warning(PRECICE_LOG_LOCATION, precice::utils::format_or_error(__VA_ARGS__))

#define PRECICE_INFO(...) _log.info(PRECICE_LOG_LOCATION, precice::utils::format_or_error(__VA_ARGS__))

#define PRECICE_ERROR(...)                                                          \
  do {                                                                              \
    _log.error(PRECICE_LOG_LOCATION, precice::utils::format_or_error(__VA_ARGS__)); \
    std::exit(-1);                                                                  \
  } while (false)

#define PRECICE_WARN_IF(condition, ...) \
  do {                                  \
    if (condition) {                    \
      PRECICE_WARN(__VA_ARGS__);        \
    }                                   \
  } while (false)

#define PRECICE_INFO_IF(condition, ...) \
  do {                                  \
    if (condition) {                    \
      PRECICE_INFO(__VA_ARGS__);        \
    }                                   \
  } while (false)

#define PRECICE_CHECK(check, ...) \
  do {                            \
    if (!(check)) {               \
      PRECICE_ERROR(__VA_ARGS__); \
    }                             \
  } while (false)

// Debug logging is disabled in release (NDEBUG) builds by default.
// To enable it anyhow, enable the CMake option PRECICE_RELEASE_WITH_DEBUG_LOG.

#if defined(NDEBUG) && !defined(PRECICE_RELEASE_WITH_DEBUG_LOG)
#define PRECICE_NO_DEBUG_LOG
#endif

#ifdef PRECICE_NO_DEBUG_LOG

#include "utils/ignore.hpp"

#define PRECICE_DEBUG(...) \
  ::precice::utils::ignore(__VA_ARGS__)

#define PRECICE_DEBUG_IF(...) \
  ::precice::utils::ignore(__VA_ARGS__)

#define PRECICE_TRACE(...) \
  ::precice::utils::ignore(__VA_ARGS__)

#else // PRECICE_NO_DEBUG_LOG

#define PRECICE_DEBUG(...) _log.debug(PRECICE_LOG_LOCATION, precice::utils::format_or_error(__VA_ARGS__))

#define PRECICE_DEBUG_IF(condition, ...) \
  do {                                   \
    if (condition) {                     \
      PRECICE_DEBUG(__VA_ARGS__);        \
    }                                    \
  } while (false)

#endif // ! PRECICE_NO_DEBUG_LOG

// Trace logging is disabled in release (NDEBUG) builds by default.
// To enable it anyhow, enable the CMake option PRECICE_RELEASE_WITH_TRACE_LOG.

#if defined(NDEBUG) && !defined(PRECICE_RELEASE_WITH_TRACE_LOG)
#define PRECICE_NO_TRACE_LOG
#endif

#ifdef PRECICE_NO_TRACE_LOG

#include "utils/ignore.hpp"

#define PRECICE_TRACE(...) \
  ::precice::utils::ignore(__VA_ARGS__)

#else // PRECICE_NO_TRACE_LOG

#include "logging/Tracer.hpp"
#include "utils/ArgumentFormatter.hpp"

// Do not put do {...} while (false) here, it will destroy the _tracer_ right after creation
#define PRECICE_TRACE(...)                                       \
  precice::logging::Tracer _tracer_(_log, PRECICE_LOG_LOCATION); \
  _log.trace(PRECICE_LOG_LOCATION, std::string{"Entering "} + __func__ + PRECICE_LOG_ARGUMENTS(__VA_ARGS__))

#endif // ! PRECICE_NO_TRACE_LOG
