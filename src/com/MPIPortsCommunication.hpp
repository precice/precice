// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_COMMUNICATIONMPIPORTS_HPP_
#define PRECICE_COM_COMMUNICATIONMPIPORTS_HPP_

#ifndef PRECICE_NO_MPI

#include "MPICommunication.hpp"
#include "tarch/logging/Log.h"

namespace precice {
namespace com {

/**
 * @brief Provides connection methods based on MPI ports (part of MPI 2.0).
 *
 * The two participants to be connected can be run in two process groups started
 * up individually, i.e. not within the same process group.
 */
class MPIPortsCommunication : public MPICommunication
{
public:

   /**
    * @brief Constructor.
    */
   MPIPortsCommunication ( const std::string& publishingDirectory );

   /**
    * @brief Destructor.
    */
   virtual ~MPIPortsCommunication();

   /**
    * @brief Returns true, if a connection to a remote participant has been setup.
    */
   virtual bool isConnected()
   {
      return _isConnection;
   }

   /**
    * @brief Returns the number of processes in the remote communicator.
    *
    * Precondition: a connection to the remote participant has been setup.
    */
   virtual int getRemoteCommunicatorSize();

   /**
    * @brief See precice::com::Communication::acceptConnection().
    */
   virtual void acceptConnection (
      const std::string& nameAcceptor,
      const std::string& nameRequester,
      int                acceptorProcessRank,
      int                acceptorCommunicatorSize );

   /**
    * @brief See precice::com::Communication::requestConnection().
    */
   virtual void requestConnection (
      const std::string& nameAcceptor,
      const std::string& nameRequester,
      int                requesterProcessRank,
      int                requesterCommunicatorSize );

   /**
    * @brief See precice::com::Communication::closeConnection().
    */
   virtual void closeConnection();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   std::string _publishingDirectory;

   // @brief Name of the port used for connection.
   char _portname[MPI_MAX_PORT_NAME];

   // @brief Flag indicating a connection.
   bool _isConnection;
};

}} // close namespaces


#endif // not PRECICE_NO_MPI
#endif /* PRECICE_COM_COMMUNICATIONMPIPORTS_HPP_ */
