#pragma once

#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/facilities/empty.hpp>

#include <boost/vmd/is_empty.hpp>

#include <string>
#include "utils/String.hpp"
#include "utils/prettyprint.hpp" // so that we can put std::vector et. al. on ostream
#include "utils/SignalHandler.hpp"

#include "Tracer.hpp"

#define WARN(message) _log.warning(LOG_LOCATION, PRECICE_AS_STRING(message))

#define INFO(message) _log.info(LOG_LOCATION, PRECICE_AS_STRING(message))

#define ERROR(message) do {                                             \
    _log.error(LOG_LOCATION, PRECICE_AS_STRING(message));               \
    precice::utils::terminationSignalHandler(6);                        \
    std::exit(-1);                                                      \
} while (false)

#define CHECK(check, message)                      \
  if ( !(check) ) {                                \
    ERROR(message);                                \
  }

#ifdef NDEBUG 

#define DEBUG(...) {}
#define TRACE(...) {}

#else // NDEBUG

#define DEBUG(message) _log.debug(LOG_LOCATION, PRECICE_AS_STRING(message))

/// Helper macro, used by TRACE
#define LOG_ARGUMENT(r, data, i, elem)                                  \
  << '\n' << "  Argument " << i << ": " << BOOST_PP_STRINGIZE(elem) << " == " << elem

// Do not put do {...} while (false) here, it will destroy the _tracer_ right after creation
#define TRACE(...)                                                      \
  /*BOOST_LOG_FUNCTION();*/                                             \
  precice::logging::Tracer _tracer_(_log, LOG_LOCATION);                \
  _log.trace(LOG_LOCATION, PRECICE_AS_STRING("Entering " << __func__    \
  BOOST_PP_IF(BOOST_VMD_IS_EMPTY(__VA_ARGS__),                          \
              BOOST_PP_EMPTY(),                                         \
              BOOST_PP_SEQ_FOR_EACH_I(LOG_ARGUMENT,, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))));


#endif // ! NDEBUG


#define LOG_LOCATION precice::logging::LogLocation{__FILE__,  __LINE__, __func__}

  







