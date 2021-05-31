#pragma once

#include "utils/fmt.hpp"

#define PRECICE_LOG_LOCATION     \
  precice::logging::LogLocation  \
  {                              \
    __FILE__, __LINE__, __func__ \
  }

#define PRECICE_WARN(...) _log.warning(PRECICE_LOG_LOCATION, fmt::format(__VA_ARGS__))

#define PRECICE_INFO(...) _log.info(PRECICE_LOG_LOCATION, fmt::format(__VA_ARGS__))

#define PRECICE_ERROR(...)                                      \
  do {                                                          \
    _log.error(PRECICE_LOG_LOCATION, fmt::format(__VA_ARGS__)); \
    std::exit(-1);                                              \
  } while (false)

#define PRECICE_CHECK(check, ...) \
  if (!(check)) {                 \
    PRECICE_ERROR(__VA_ARGS__);   \
  }

#ifdef NDEBUG

#define PRECICE_DEBUG(...) \
  {                        \
  }
#define PRECICE_TRACE(...) \
  {                        \
  }

#else // NDEBUG

#include "logging/Tracer.hpp"
#include "utils/ArgumentFormatter.hpp"

#define PRECICE_DEBUG(...) _log.debug(PRECICE_LOG_LOCATION, fmt::format(__VA_ARGS__))

// Do not put do {...} while (false) here, it will destroy the _tracer_ right after creation
#define PRECICE_TRACE(...)                                       \
  precice::logging::Tracer _tracer_(_log, PRECICE_LOG_LOCATION); \
  _log.trace(PRECICE_LOG_LOCATION, std::string{"Entering "} + __func__ + PRECICE_LOG_ARGUMENTS(__VA_ARGS__))

#endif // ! NDEBUG
