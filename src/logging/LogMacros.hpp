#pragma once

#include <boost/log/expressions.hpp>
#include <boost/log/attributes/mutable_constant.hpp>

#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/facilities/empty.hpp>

#include <boost/vmd/is_empty.hpp>

#include <string>
  
#define WARN(message) do {                                  \
    LOG_LOCATION;                                                       \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::warning)   \
      << message;                                                       \
  } while (false)

#define INFO(message)                                                   \
  if (not precice::utils::MasterSlave::_slaveMode) {                    \
    LOG_LOCATION;                                                       \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::info)      \
      << message;                                                       \
  }

#define ERROR(message) do {                                             \
    LOG_LOCATION;                                                       \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::error)     \
      << message;                                                       \
    std::abort();                                                       \
  } while (false)

#define CHECK(check, message)          \
  if ( !(check) ) {                                \
    ERROR(message);                                \
  }

#define preciceInfo(methodname, message) INFO(message)
#define preciceWarning(methodname, message) WARN(message)
#define preciceError(methodname, message) ERROR(message)
#define preciceCheck(check, methodname, message) CHECK(check, message)

#ifdef NDEBUG 

#define DEBUG(...)
#define TRACE(...)
#define preciceTrace(...)

#else // NDEBUG

#define DEBUG(message) do {                                             \
    LOG_LOCATION;                                                       \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::debug)     \
      << message;                                                       \
  } while (false)

#include "Tracer.hpp"

/// Helper macro, used by preciceTrace.
#define LOG_ARGUMENT(r, data, i, elem)                  \
  << std::endl << "  Argument " << i << ": " << BOOST_PP_STRINGIZE(elem) << " == " << elem


// Do not put do {...} while (false) here, it will destroy the _tracer_ right after creation
#define TRACE(...)                                                      \
  LOG_LOCATION;                                                         \
  BOOST_LOG_FUNCTION();                                                 \
  precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__); \
  BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace) << "Entering " << __func__ \
  BOOST_PP_IF(BOOST_VMD_IS_EMPTY(__VA_ARGS__),                          \
              BOOST_PP_EMPTY(),                                         \
              BOOST_PP_SEQ_FOR_EACH_I(LOG_ARGUMENT,, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)));

/// Legacy macro, removes the first element of __VA_ARGS__
#define preciceTrace(...)                                                      \
  LOG_LOCATION;                                                         \
  BOOST_LOG_FUNCTION();                                                 \
  precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__); \
  BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace) << "Entering " << __func__ \
  BOOST_PP_SEQ_FOR_EACH_I(LOG_ARGUMENT,, BOOST_PP_SEQ_TAIL(BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)));
  


#endif // ! NDEBUG


#define LOG_LOCATION do {                                               \
    boost::log::attribute_cast<boost::log::attributes::mutable_constant<int>>( \
      boost::log::core::get()->get_global_attributes()["Line"]).set(__LINE__); \
    boost::log::attribute_cast<boost::log::attributes::mutable_constant<std::string>>( \
      boost::log::core::get()->get_global_attributes()["File"]).set(__FILE__); \
    boost::log::attribute_cast<boost::log::attributes::mutable_constant<std::string>>( \
      boost::log::core::get()->get_global_attributes()["Function"]).set(__func__); \
  } while (false)

  







