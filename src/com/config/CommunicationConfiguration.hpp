// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_COMMUNICATIONCONFIGURATION_HPP_
#define PRECICE_COM_COMMUNICATIONCONFIGURATION_HPP_

#include "com/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "utils/xml/XMLTag.hpp"
#include <string>
#include <vector>
#include "boost/tuple/tuple.hpp"

namespace precice {
namespace com {

/**
 * @brief Configuration for communication channels between solvers.
 */
class CommunicationConfiguration : public utils::XMLTag::Listener
{
public:

   typedef boost::tuple<PtrCommunication,std::string,std::string> ComTuple;

   /**
    * @brief Returns the name of this configurations enclosing tag.
    */
   //static const std::string& getTag();

   /**
    * @brief Creates a not auto-configurable config, to use createCommunicatio().
    */
   CommunicationConfiguration();

   /**
    * @brief Constructor.
    */
   CommunicationConfiguration ( utils::XMLTag& parent );

   //bool parseSubtag ( utils::XMLTag::XMLReader* xmlReader );

//   bool isValid() const
//   {
//      return _isValid;
//   }

   /**
    * @brief Returns the communication object for the given user names.
    *
    * Exits with an error message, when no object is configured for the given
    * user names.
    */
   PtrCommunication getCommunication (
      const std::string& from,
      const std::string& to );

   /**
    * @brief Returns all configured communication objects.
    */
   std::vector<ComTuple>& communications()
   {
      return _communications;
   }

   /**
    * @brief Adds a communication object of given type. For manual configuration.
    */
//   void addCommunication (
//      const std::string& type,
//      const std::string& from,
//      const std::string& to,
//      const std::string& context );

   virtual void xmlTagCallback ( utils::XMLTag& callingTag );

   virtual void xmlEndTagCallback ( utils::XMLTag& callingTag ) {}

   /**
    * @brief Returns a communication object of given type.
    */
   PtrCommunication createCommunication ( const utils::XMLTag& tag ) const;

private:

   static tarch::logging::Log _log;

   const std::string TAG;
   const std::string ATTR_TYPE;
   const std::string ATTR_FROM;
   const std::string ATTR_TO;
   const std::string ATTR_PORT;
   const std::string ATTR_EXCHANGE_DIRECTORY;

   const std::string VALUE_MPI;
   const std::string VALUE_MPI_SINGLE;
   const std::string VALUE_FILES;
   const std::string VALUE_SOCKETS;

   std::vector<ComTuple> _communications;

   //bool _isValid;

   void checkDuplicates (
     const std::string& from,
     const std::string& to );
};

}} // namespace precice, com

#endif /* PRECICE_COM_COMMUNICATIONCONFIGURATION_HPP_ */
