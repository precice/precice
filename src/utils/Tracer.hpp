// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_TRACER_HPP_
#define PRECICE_UTILS_TRACER_HPP_

#include "tarch/logging/Log.h"
#include <string>

namespace precice {
namespace utils {

/**
 * @brief Logging assistant for trace logging into and out of methods.
 *
 * The tracer creates an entry log on creation and an exit log on destruction,
 * i.e. when the scope of the creation is left. In combination with the macro
 * preciceTrace(), defined in utils/Helpers.hpp, a convenient tracing mechanism
 * is provided.
 */
class Tracer
{
public:

   /**
    * @brief Constructor, logs entry message.
    *
    * @param log [IN] Used to create the entry and exit logs.
    * @param methodname [IN] Name of method to be traced.
    * @param stateString [IN] Additional string printed with entry log.
    */
   Tracer (
     const tarch::logging::Log & log,
     const std::string & methodname,
     const std::string & stateString );

   /**
    * @brief Destructor, logs exit message.
    */
   ~Tracer();

private:

   // @brief Logging device, called _log to enable use of logging macros.
   const tarch::logging::Log & _log;

   // @brief Name of the method to be traced.
   std::string _methodname;

   // @brief State string printed with entry log.
   std::string _stateString;
};

}} // namespace precice, utils

#endif /* PRECICE_UTILS_TRACER_HPP_ */
