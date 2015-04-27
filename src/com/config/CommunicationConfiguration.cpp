// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CommunicationConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/FileCommunication.hpp"
#include "m2n/M2N.hpp"
#ifndef PRECICE_NO_SOCKETS
#include "com/SocketCommunication.hpp"
# endif // not PRECICE_NO_SOCKETS
#include "utils/Globals.hpp"
#include "utils/xml/XMLAttribute.hpp"

namespace precice {
namespace com {

tarch::logging::Log CommunicationConfiguration::
   _log("precice::com::CommunicationConfiguration");

CommunicationConfiguration:: CommunicationConfiguration()
:
  TAG("communication"),
  ATTR_TYPE("type"),
  ATTR_FROM("from"),
  ATTR_TO("to"),
  ATTR_PORT("port"),
  ATTR_NETWORK("network"),
  ATTR_EXCHANGE_DIRECTORY("exchange-directory"),
  VALUE_MPI("mpi"),
  VALUE_MPI_SINGLE("mpi-single"),
  VALUE_FILES("files"),
  VALUE_SOCKETS("sockets")
{}

Communication::SharedPointer CommunicationConfiguration:: createCommunication
(
  const utils::XMLTag& tag ) const
{
  com::Communication::SharedPointer com;
  if (tag.getName() == VALUE_SOCKETS){
#   ifdef PRECICE_NO_SOCKETS
    std::ostringstream error;
    error << "Communication type \"" << VALUE_SOCKETS << "\" can only be used "
          << "when preCICE is compiled with argument \"sockets=on\"";
    throw error.str();
#   else
    std::string network = tag.getStringAttributeValue(ATTR_NETWORK);
    int port = tag.getIntAttributeValue(ATTR_PORT);

    preciceCheck(not utils::isTruncated<unsigned short>(port),
                 "createCommunication()",
                 "The value given for the \"port\" attribute is not a "
                 "16-bit unsigned integer: " << port);

    std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
    com = com::Communication::SharedPointer(new com::SocketCommunication(port, false, network, dir));
#   endif // PRECICE_NO_SOCKETS
  }
  else if (tag.getName() == VALUE_MPI){
    std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
#   ifdef PRECICE_NO_MPI
    std::ostringstream error;
    error << "Communication type \"" << VALUE_MPI << "\" can only be used "
          << "when preCICE is compiled with argument \"mpi=on\"";
    throw error.str();
#   else
    com = com::Communication::SharedPointer(new com::MPIPortsCommunication(dir));
#   endif
  }
  else if (tag.getName() == VALUE_MPI_SINGLE){
#   ifdef PRECICE_NO_MPI
    std::ostringstream error;
    error << "Communication type \"" << VALUE_MPI_SINGLE << "\" can only be used "
          << "when preCICE is compiled with argument \"mpi=on\"";
    throw error.str();
#   else
    com = com::Communication::SharedPointer(new com::MPIDirectCommunication());
#   endif
  }
  else if (tag.getName() == VALUE_FILES){
    std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
    com = com::Communication::SharedPointer(new com::FileCommunication(false, dir));
  }
  assertion(com.get() != NULL);
  return com;
}

}} // namespace precice, com
