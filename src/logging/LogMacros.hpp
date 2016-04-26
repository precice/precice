#pragma once

//#include "../utils/MasterSlave.hpp"

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/expressions.hpp>
#include <boost/log/attributes/mutable_constant.hpp>

#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

  
#define precicePrint(message) do                                            \
  {                                                                          \
    preciceInfo("unknown",message);                                         \
  } while (false)


#define preciceDebug(message) do                                           \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::debug)         \
    << message;                                                             \
  } while (false)

    
#define preciceWarning(methodname, message) do                             \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::warning)       \
    << message;                                                             \
  } while (false)
    
#define preciceInfo(methodname, message) do                                \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::info)          \
    << message;                                                             \
  } while (false)

#ifdef NDEBUG 

//#define preciceDebug(methodname, message)
#define preciceTrace(...)

#else // NDEBUG

#include "Tracer.hpp"

/// Helper macro, used by preciceTrace.
#define LOG_ARGUMENT(r, data, i, elem)                        \
  << "  Argument " << i << ": " << elem << std::endl

#define preciceTrace(...) do { \
    LOG_LOCATION;                                                       \
    BOOST_LOG_FUNCTION();                                               \
    precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__); \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace) << "Entering " << __func__ \
    BOOST_PP_SEQ_FOR_EACH_I(LOG_ARGUMENT,, BOOST_PP_SEQ_TAIL(BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))); \
  } while (false)

#endif // ! NDEBUG

#define preciceTrace1 preciceTrace
#define preciceTrace2 preciceTrace
#define preciceTrace3 preciceTrace
#define preciceTrace4 preciceTrace
#define preciceTrace5 preciceTrace
#define preciceTrace6 preciceTrace

#define preciceError(methodname, message) do                             \
  {                                                                       \
    LOG_LOCATION;                                                         \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::error)       \
    << message;                                                            \
    std::abort();                                                         \
  } while (false)

#define preciceCheck(check, methodname, errormessage)                    \
  if ( !(check) ) {                                                       \
    preciceError(methodname, errormessage);                              \
  }

  
#define LOG_LOCATION do                                                                 \
  {                                                                                     \
    boost::log::attribute_cast<boost::log::attributes::mutable_constant<int>>(          \
      boost::log::core::get()->get_global_attributes()["Line"]).set(__LINE__);          \
    boost::log::attribute_cast<boost::log::attributes::mutable_constant<std::string>>(  \
      boost::log::core::get()->get_global_attributes()["File"]).set(__FILE__);          \
    boost::log::attribute_cast<boost::log::attributes::mutable_constant<std::string>>(  \
      boost::log::core::get()->get_global_attributes()["Function"]).set(__func__);      \
  } while (false)

  







