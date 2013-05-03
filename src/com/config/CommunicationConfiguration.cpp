// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CommunicationConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/FileCommunication.hpp"
#ifndef PRECICE_NO_SOCKETS
#include "com/SocketCommunication.hpp"
# endif // not PRECICE_NO_SOCKETS
#include "utils/Globals.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"

namespace precice {
namespace com {

tarch::logging::Log CommunicationConfiguration::
   _log("precice::com::CommunicationConfiguration");

//const std::string & CommunicationConfiguration:: getTag()
//{
//  static std::string tag("communication");
//  return tag;
//}

CommunicationConfiguration:: CommunicationConfiguration()
:
  TAG("communication"),
  ATTR_TYPE("type"),
  ATTR_FROM("from"),
  ATTR_TO("to"),
  ATTR_PORT("port"),
  ATTR_EXCHANGE_DIRECTORY("exchange-directory"),
  VALUE_MPI("mpi"),
  VALUE_MPI_SINGLE("mpi-single"),
  VALUE_FILES("files"),
  VALUE_SOCKETS("sockets"),
  _communications()
  //_isValid(false)
{}

CommunicationConfiguration:: CommunicationConfiguration
(
  utils::XMLTag& parent )
:
  TAG("communication"),
  ATTR_TYPE("type"),
  ATTR_FROM("from"),
  ATTR_TO("to"),
  ATTR_PORT("port"),
  ATTR_EXCHANGE_DIRECTORY("exchange-directory"),
  VALUE_MPI("mpi"),
  VALUE_MPI_SINGLE("mpi-single"),
  VALUE_FILES("files"),
  VALUE_SOCKETS("sockets"),
  _communications()
  //_isValid(false)
{
  using namespace utils;
  std::string doc;
  std::list<XMLTag> tags;
  XMLTag::Occurrence occ = XMLTag::OCCUR_NOT_OR_ONCE;
  {
    XMLTag tag(*this, VALUE_SOCKETS, occ, TAG);
    doc = "Communication via Sockets.";
    tag.setDocumentation(doc);

    XMLAttribute<int> attrPort(ATTR_PORT);
    doc = "Port number to be used by server for socket communiation.";
    attrPort.setDocumentation(doc);
    attrPort.setDefaultValue(51235);
    tag.addAttribute(attrPort);

    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_MPI, occ, TAG);
    doc = "Communication via MPI with startup in separated communication spaces.";
    tag.setDocumentation(doc);

    XMLAttribute<std::string> attrExchangeDirectory(ATTR_EXCHANGE_DIRECTORY);
    doc = "Directory where connection information is exchanged. By default, the ";
    doc += "directory of startup is chosen, and both solvers have to be started ";
    doc += "in the same directory.";
    attrExchangeDirectory.setDocumentation(doc);
    attrExchangeDirectory.setDefaultValue("");
    tag.addAttribute(attrExchangeDirectory);

    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_MPI_SINGLE, occ, TAG);
    doc = "Communication via MPI with startup in common communication space.";
    tag.setDocumentation(doc);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_FILES, occ, TAG);
    doc = "Communication via files.";
    tag.setDocumentation(doc);

    XMLAttribute<std::string> attrExchangeDirectory(ATTR_EXCHANGE_DIRECTORY);
    doc = "Directory where communication files are exchanged. By default, the ";
    doc += "directory of startup is chosen, and both solvers have to be started ";
    doc += "in the same directory.";
    attrExchangeDirectory.setDocumentation(doc);
    attrExchangeDirectory.setDefaultValue("");
    tag.addAttribute(attrExchangeDirectory);

    tags.push_back(tag);
  }

//  XMLAttribute<std::string> attrType(ATTR_TYPE);
//  ValidatorEquals<std::string> validMPI(VALUE_MPI);
//  ValidatorEquals<std::string> validMPISingle(VALUE_MPI_SINGLE);
//  ValidatorEquals<std::string> validFiles(VALUE_FILES);
//  ValidatorEquals<std::string> validSockets(VALUE_SOCKETS);
//  attrType.setValidator(validMPI || validMPISingle || validFiles || validSockets);
//  tag.addAttribute(attrType);

  XMLAttribute<std::string> attrFrom ( ATTR_FROM );
  doc = "First participant name involved in communication.";
  attrFrom.setDocumentation(doc);
  XMLAttribute<std::string> attrTo(ATTR_TO);
  doc = "Second participant name involved in communication.";
  attrTo.setDocumentation(doc);

  foreach (XMLTag& tag, tags){
    tag.addAttribute(attrFrom);
    tag.addAttribute(attrTo);
    parent.addSubtag(tag);
  }
}

//bool CommunicationConfiguration:: parseSubtag
//(
//  utils::XMLTag::XMLReader* xmlReader )
//{
//  preciceTrace ( "parseSubtag()" );
//  using namespace utils;
//  XMLTag tag ( TAG, XMLTag::OCCUR_ONCE );
//
//  XMLAttribute<std::string> attrType ( ATTR_TYPE );
//  ValidatorEquals<std::string> validMPI ( VALUE_MPI );
//  ValidatorEquals<std::string> validMPISingle ( VALUE_MPI_SINGLE );
//  ValidatorEquals<std::string> validFiles ( VALUE_FILES );
//  ValidatorEquals<std::string> validSockets ( VALUE_SOCKETS );
//  attrType.setValidator ( validMPI || validMPISingle || validFiles || validSockets );
//  tag.addAttribute(attrType);
//
//  XMLAttribute<std::string> attrFrom ( ATTR_FROM );
//  tag.addAttribute(attrFrom);
//
//  XMLAttribute<std::string> attrTo ( ATTR_TO );
//  tag.addAttribute(attrTo);
//
//  XMLAttribute<std::string> attrContext(ATTR_CONTEXT);
//  attrContext.setDefaultValue("");
//  tag.addAttribute(attrContext);
//
//  _isValid = tag.parse(xmlReader, *this);
//  return _isValid;
//}

PtrCommunication CommunicationConfiguration:: getCommunication
(
  const std::string& from,
  const std::string& to )
{
  using boost::get;
  foreach (ComTuple & tuple, _communications){
    if ((get<1>(tuple) == from) && (get<2>(tuple) == to)){
      return get<0>(tuple);
    }
    else if ((get<2>(tuple) == from) && (get<1>(tuple) == to)){
      return get<0>(tuple);
    }
  }
  std::ostringstream error;
  error << "No communication configured between \"" << from << "\" and \""
        << to << "\"!";
  throw error.str();
}

//void CommunicationConfiguration:: addCommunication
//(
//  const std::string& type,
//  const std::string& from,
//  const std::string& to,
//  const std::string& context )
//{
//  using boost::get;
//  bool alreadyAdded = false;
//  foreach (ComTuple& tuple, _communications){
//    alreadyAdded |= (get<1>(tuple) == from) && (get<2>(tuple) == to);
//    alreadyAdded |= (get<2>(tuple) == from) && (get<1>(tuple) == to);
//  }
//  if (alreadyAdded){
//    std::ostringstream stream;
//    stream << "Multiple communication between between user \"" << from
//           << "\" and \"" << to << "\" is invalid!";
//    throw stream.str();
//  }
//  _communications.push_back (
//      boost::make_tuple(createCommunication(type, context), from, to) );
//}

void CommunicationConfiguration:: xmlTagCallback
(
   utils::XMLTag& tag )
{
  if (tag.getNamespace() == TAG){
    std::string from = tag.getStringAttributeValue(ATTR_FROM);
    std::string to = tag.getStringAttributeValue(ATTR_TO);
    checkDuplicates(from, to);
    com::PtrCommunication com = createCommunication(tag);
    assertion(com.get() != NULL);
    _communications.push_back(boost::make_tuple(com, from, to));
  }
}

PtrCommunication CommunicationConfiguration:: createCommunication
(
  const utils::XMLTag& tag ) const
{
  com::PtrCommunication com;
  //std::string from = tag.getStringAttributeValue(ATTR_FROM);
  //std::string to = tag.getStringAttributeValue(ATTR_TO);
  if (tag.getName() == VALUE_SOCKETS){
#   ifdef PRECICE_NO_SOCKETS
    std::ostringstream error;
    error << "Communication type \"" << VALUE_SOCKETS << "\" can only be used "
          << "when preCICE is compiled with argument \"sockets=on\"";
    throw error.str();
#   else
    int port = tag.getIntAttributeValue(ATTR_PORT);
    com = com::PtrCommunication(new com::SocketCommunication(port));
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
    com = com::PtrCommunication(new com::MPIPortsCommunication(dir));
#   endif
  }
  else if (tag.getName() == VALUE_MPI_SINGLE){
#   ifdef PRECICE_NO_MPI
    std::ostringstream error;
    error << "Communication type \"" << VALUE_MPI_SINGLE << "\" can only be used "
          << "when preCICE is compiled with argument \"mpi=on\"";
    throw error.str();
#   else
    com = com::PtrCommunication(new com::MPIDirectCommunication());
#   endif
  }
  else if (tag.getName() == VALUE_FILES){
    std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
    com = com::PtrCommunication(new com::FileCommunication(false, dir));
  }
  assertion(com.get() != NULL);
  return com;
}

void CommunicationConfiguration:: checkDuplicates
(
  const std::string& from,
  const std::string& to )
{
  using boost::get;
  bool alreadyAdded = false;
  foreach (ComTuple& tuple, _communications){
    alreadyAdded |= (get<1>(tuple) == from) && (get<2>(tuple) == to);
    alreadyAdded |= (get<2>(tuple) == from) && (get<1>(tuple) == to);
  }
  if (alreadyAdded){
    std::ostringstream error;
    error << "Multiple communication defined between participant \"" << from
          << "\" and \"" << to << "\"";
    throw error.str();
  }
}

}} // namespace precice, com
