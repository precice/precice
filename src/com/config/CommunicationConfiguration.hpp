#pragma once

#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "utils/xml/XMLTag.hpp"

#include <string>

namespace precice {
namespace com {

/**
 * @brief Configuration for communication channels between server and clients or master and slaves.
 * The communication between two solvers is configured in m2n::M2NConfiguration
 */
class CommunicationConfiguration
{
public:

   /**
    * @brief Constructor
    */
   CommunicationConfiguration();

   virtual ~CommunicationConfiguration() {}

   /**
    * @brief Returns a communication object of given type.
    */
   PtrCommunication createCommunication ( const utils::XMLTag& tag ) const;

private:

   static logging::Logger _log;

   const std::string TAG;
   const std::string ATTR_TYPE;
   const std::string ATTR_FROM;
   const std::string ATTR_TO;
   const std::string ATTR_PORT;
   const std::string ATTR_NETWORK;
   const std::string ATTR_EXCHANGE_DIRECTORY;

   const std::string VALUE_MPI;
   const std::string VALUE_MPI_SINGLE;
   const std::string VALUE_SOCKETS;

};

}} // namespace precice, com
