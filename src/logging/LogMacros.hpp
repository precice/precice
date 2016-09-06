#pragma once

#include <boost/log/expressions.hpp>
#include <boost/log/attributes/mutable_constant.hpp>

#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>
  
#define precicePrint(message) do {              \
    preciceInfo("unknown", message);            \
  } while (false)

#define preciceWarning(methodname, message) do {                        \
    LOG_LOCATION;                                                       \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::warning)   \
      << message;                                                       \
  } while (false)
    
#define preciceInfo(methodname, message)                                \
  if (not precice::utils::MasterSlave::_slaveMode) {                      \
    LOG_LOCATION;                                                       \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::info)      \
      << message;                                                       \
  } 

#ifdef NDEBUG 

//#define preciceDebug(methodname, message)
#define preciceDebug(...)
#define preciceTrace(...)

#else // NDEBUG

#include "Tracer.hpp"

/// Helper macro, used by preciceTrace.
#define LOG_ARGUMENT(r, data, i, elem)                  \
  << std::endl << "  Argument " << i << ": " << elem

#define preciceDebug(message) do {                                      \
    LOG_LOCATION;                                                       \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::debug)     \
      << message;                                                       \
  } while (false)

// Do not put do {...} while (false) here, it will destroy the _tracer_ right after creation
#define preciceTrace(...)                                               \
  LOG_LOCATION;                                                         \
  BOOST_LOG_FUNCTION();                                                 \
  precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__); \
  BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace) << "Entering " << __func__ \
  BOOST_PP_SEQ_FOR_EACH_I(LOG_ARGUMENT,, BOOST_PP_SEQ_TAIL(BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)));
  
#endif // ! NDEBUG

#define preciceError(methodname, message) do {                          \
    LOG_LOCATION;                                                       \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::error)     \
      << message;                                                       \
    std::abort();                                                       \
  } while (false)

#define preciceCheck(check, methodname, errormessage)   \
  if ( !(check) ) {                                     \
    preciceError(methodname, errormessage);             \
  }

  
#define LOG_LOCATION do {                                               \
    boost::log::attribute_cast<boost::log::attributes::mutable_constant<int>>( \
      boost::log::core::get()->get_global_attributes()["Line"]).set(__LINE__); \
    boost::log::attribute_cast<boost::log::attributes::mutable_constant<std::string>>( \
      boost::log::core::get()->get_global_attributes()["File"]).set(__FILE__); \
    boost::log::attribute_cast<boost::log::attributes::mutable_constant<std::string>>( \
      boost::log::core::get()->get_global_attributes()["Function"]).set(__func__); \
  } while (false)

  







