#pragma once

//#include "../utils/MasterSlave.hpp"

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/expressions.hpp>
#include <boost/log/attributes/mutable_constant.hpp>

#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

//#define assertion assertion

  
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

#ifdef Debug

#include "Tracer.hpp"

#define preciceTrace(methodname) do                                        \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_FUNCTION();                                                   \
    precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__);   \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace)         \
    << "Entering " << __func__;                                             \
  } while (false)

#define preciceTrace1(methodname, var1) do                                 \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_FUNCTION();                                                   \
    precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__);   \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace)         \
    << "Entering " << __func__                                              \
    << "\n" << #var1 << " = " << var1;                                      \
  } while (false)

#define preciceTrace2(methodname, var1, var2) do                           \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_FUNCTION();                                                   \
    precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__);   \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace)         \
    << "Entering " << __func__                                              \
    << "\n" << #var1 << " = " << var1                                       \
    << "\n" << #var2 << " = " << var2;                                      \
  } while (false)
 
#define preciceTrace3(methodname, var1, var2, var3) do                     \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_FUNCTION();                                                   \
    precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__);   \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace)         \
    << "Entering " << __func__                                              \
    << "\n" << #var1 << " = " << var1                                       \
    << "\n" << #var2 << " = " << var2                                       \
    << "\n" << #var3 << " = " << var3;                                      \
  } while (false)

#define preciceTrace4(methodname, var1, var2, var3, var4) do               \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_FUNCTION();                                                   \
    precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__);   \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace)         \
    << "Entering " << __func__                                              \
    << "\n" << #var1 << " = " << var1                                       \
    << "\n" << #var2 << " = " << var2                                       \
    << "\n" << #var3 << " = " << var3                                       \
    << "\n" << #var4 << " = " << var4;                                      \
  } while (false)
  
#define preciceTrace5(methodname, var1, var2, var3, var4, var5) do         \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_FUNCTION();                                                   \
    precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__);   \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace)         \
    << "Entering " << __func__                                              \
    << "\n" << #var1 << " = " << var1                                       \
    << "\n" << #var2 << " = " << var2                                       \
    << "\n" << #var3 << " = " << var3                                       \
    << "\n" << #var4 << " = " << var4                                       \
    << "\n" << #var5 << " = " << var5;                                      \
  } while (false)
  
#define preciceTrace6(methodname, var1, var2, var3, var4, var5, var6) do   \
  {                                                                         \
    LOG_LOCATION;                                                           \
    BOOST_LOG_FUNCTION();                                                   \
    precice::logging::Tracer _tracer_(_log, __func__, __FILE__,__LINE__);   \
    BOOST_LOG_SEV(_log, boost::log::trivial::severity_level::trace)         \
    << "Entering " << __func__                                              \
    << "\n" << #var1 << " = " << var1                                       \
    << "\n" << #var2 << " = " << var2                                       \
    << "\n" << #var3 << " = " << var3                                       \
    << "\n" << #var4 << " = " << var4                                       \
    << "\n" << #var5 << " = " << var5                                       \
    << "\n" << #var6 << " = " << var6;                                      \
  } while (false)
  
#else // Debug

//#define preciceDebug(methodname, message)
#define preciceTrace(methodname)
#define preciceTrace1(methodname, var1)
#define preciceTrace2(methodname, var1, var2)
#define preciceTrace3(methodname, var1, var2, var3)
#define preciceTrace4(methodname, var1, var2, var3, var4)
#define preciceTrace5(methodname, var1, var2, var3, var4, var5)
#define preciceTrace6(methodname, var1, var2, var3, var4, var5, var6)

#endif // ! Debug


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

  







