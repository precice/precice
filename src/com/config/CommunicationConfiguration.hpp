#pragma once

#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace com {

/**
 * @brief Configuration for communication channels between server and clients or master and slaves.
 * The communication between two solvers is configured in m2n::M2NConfiguration
 */
class CommunicationConfiguration
{
public:

   virtual ~CommunicationConfiguration() {}

  /// Returns a communication object of given type.
  PtrCommunication createCommunication ( const xml::XMLTag& tag ) const;

private:
  
   static logging::Logger _log;      
};

}} // namespace precice, com
