#pragma once

#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/vmd/is_empty.hpp>
#include <string>

#include "logging/Tracer.hpp"
#include "prettyprint/prettyprint.hpp" // so that we can put std::vector et. al. on ostream
#include "utils/String.hpp"

#define PRECICE_WARN(message) _log.warning(PRECICE_LOG_LOCATION, PRECICE_AS_STRING(message))

#define PRECICE_INFO(message) _log.info(PRECICE_LOG_LOCATION, PRECICE_AS_STRING(message))

#define PRECICE_ERROR(message)                                    \
  do {                                                            \
    _log.error(PRECICE_LOG_LOCATION, PRECICE_AS_STRING(message)); \
    std::exit(-1);                                                \
  } while (false)

#define PRECICE_CHECK(check, message) \
  if (!(check)) {                     \
    PRECICE_ERROR(message);           \
  }

#ifdef NDEBUG

#define PRECICE_DEBUG(...) \
  {                        \
  }
#define PRECICE_TRACE(...) \
  {                        \
  }

#else // NDEBUG

#define PRECICE_DEBUG(message) _log.debug(PRECICE_LOG_LOCATION, PRECICE_AS_STRING(message))

/// Helper macro, used by TRACE
#define PRECICE_LOG_ARGUMENT(r, data, i, elem) \
  << '\n'                                      \
  << "  Argument " << i << ": " << BOOST_PP_STRINGIZE(elem) << " == " << elem

// Do not put do {...} while (false) here, it will destroy the _tracer_ right after creation
#define PRECICE_TRACE(...)                                                                                                \
  precice::logging::Tracer _tracer_(_log, PRECICE_LOG_LOCATION);                                                          \
  _log.trace(PRECICE_LOG_LOCATION, PRECICE_AS_STRING("Entering " << __func__ BOOST_PP_IF(BOOST_VMD_IS_EMPTY(__VA_ARGS__), \
                                                                                         BOOST_PP_EMPTY(),                \
                                                                                         BOOST_PP_SEQ_FOR_EACH_I(PRECICE_LOG_ARGUMENT, , BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))));

#endif // ! NDEBUG

#define PRECICE_LOG_LOCATION     \
  precice::logging::LogLocation  \
  {                              \
    __FILE__, __LINE__, __func__ \
  }
