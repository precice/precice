// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "M2NConfiguration.hpp"
#include "m2n/M2N.hpp"
#include "m2n/DistributedComFactory.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/PointToPointComFactory.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/FileCommunication.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "utils/Globals.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"

namespace precice {
namespace m2n {

tarch::logging::Log M2NConfiguration::
   _log("precice::m2n::M2NConfiguration");

M2NConfiguration:: M2NConfiguration
(
  utils::XMLTag& parent )
:
  TAG("m2n"),
  ATTR_TYPE("type"),
  ATTR_DISTRIBUTION_TYPE("distribution-type"),
  ATTR_FROM("from"),
  ATTR_TO("to"),
  ATTR_PORT("port"),
  ATTR_NETWORK("network"),
  ATTR_EXCHANGE_DIRECTORY("exchange-directory"),
  VALUE_MPI("mpi"),
  VALUE_MPI_SINGLE("mpi-single"),
  VALUE_FILES("files"),
  VALUE_SOCKETS("sockets"),
  VALUE_GATHER_SCATTER("gather-scatter"),
  VALUE_POINT_TO_POINT("point-to-point"),
  _m2ns()
{
  using namespace utils;
  std::string doc;
  std::list<XMLTag> tags;
  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  {
    XMLTag tag(*this, VALUE_SOCKETS, occ, TAG);
    doc = "Communication via Sockets.";
    tag.setDocumentation(doc);

    XMLAttribute<int> attrPort(ATTR_PORT);
    doc = "Port number (16-bit unsigned integer) to be used for socket ";
    doc += "communiation. The default is \"0\", what means that OS will ";
    doc += "dynamically search for a free port (if at least one exists) and ";
    doc += "bind it automatically.";
    attrPort.setDocumentation(doc);
    attrPort.setDefaultValue(0);
    tag.addAttribute(attrPort);

    XMLAttribute<std::string> attrNetwork(ATTR_NETWORK);
    doc = "Network name to be used for socket communiation. ";
    doc += "Default is \"lo\", i.e., the local host loopback.";
    attrNetwork.setDocumentation(doc);
    attrNetwork.setDefaultValue("lo");
    tag.addAttribute(attrNetwork);

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


  XMLAttribute<std::string> attrDistrTypeBoth ( ATTR_DISTRIBUTION_TYPE);
  doc = "Distribution manner of the M2N communication .";
  attrDistrTypeBoth.setDocumentation(doc);
  ValidatorEquals<std::string> validDistrGatherScatter ( VALUE_GATHER_SCATTER );
  ValidatorEquals<std::string> validDistrP2P ( VALUE_POINT_TO_POINT);
  attrDistrTypeBoth.setValidator ( validDistrGatherScatter || validDistrP2P );
  attrDistrTypeBoth.setDefaultValue(VALUE_POINT_TO_POINT);

  XMLAttribute<std::string> attrDistrTypeOnly ( ATTR_DISTRIBUTION_TYPE);
  doc = "Distribution manner of the M2N communication .";
  attrDistrTypeOnly.setDocumentation(doc);
  attrDistrTypeOnly.setValidator ( validDistrGatherScatter );
  attrDistrTypeOnly.setDefaultValue(VALUE_GATHER_SCATTER);

  XMLAttribute<std::string> attrFrom ( ATTR_FROM );
  doc = "First participant name involved in communication.";
  attrFrom.setDocumentation(doc);
  XMLAttribute<std::string> attrTo(ATTR_TO);
  doc = "Second participant name involved in communication.";
  attrTo.setDocumentation(doc);

  foreach (XMLTag& tag, tags){
    tag.addAttribute(attrFrom);
    tag.addAttribute(attrTo);
    if(tag.getName() == VALUE_MPI || tag.getName() == VALUE_SOCKETS){
      tag.addAttribute(attrDistrTypeBoth);
    }
    else{
      tag.addAttribute(attrDistrTypeOnly);
    }
    parent.addSubtag(tag);
  }
}

m2n::M2N::SharedPointer M2NConfiguration:: getM2N
(
  const std::string& from,
  const std::string& to )
{
  using boost::get;
  for (M2NTuple & tuple : _m2ns){
    if ((get<1>(tuple) == from) && (get<2>(tuple) == to)){
      return get<0>(tuple);
    }
    else if ((get<2>(tuple) == from) && (get<1>(tuple) == to)){
      return get<0>(tuple);
    }
  }
  std::ostringstream error;
  error << "No m2n communication configured between \"" << from << "\" and \""
        << to << "\"!";
  throw error.str();
}

void M2NConfiguration:: xmlTagCallback
(
   utils::XMLTag& tag )
{
  if (tag.getNamespace() == TAG){
    std::string from = tag.getStringAttributeValue(ATTR_FROM);
    std::string to = tag.getStringAttributeValue(ATTR_TO);
    checkDuplicates(from, to);
    std::string distrType = tag.getStringAttributeValue(ATTR_DISTRIBUTION_TYPE);

    com::CommunicationFactory::SharedPointer comFactory;
    com::Communication::SharedPointer com;
    if (tag.getName() == VALUE_SOCKETS){
#     ifdef PRECICE_NO_SOCKETS
        std::ostringstream error;
        error << "M2N type \"" << VALUE_SOCKETS << "\" can only be used "
              << "when preCICE is compiled with argument \"sockets=on\"";
        throw error.str();
#     else
        std::string network = tag.getStringAttributeValue(ATTR_NETWORK);
        int port = tag.getIntAttributeValue(ATTR_PORT);

        preciceCheck(not utils::isTruncated<unsigned short>(port),
                     "xmlTagCallback()",
                     "The value given for the \"port\" attribute is not a "
                     "16-bit unsigned integer: " << port);

        std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
        comFactory = com::CommunicationFactory::SharedPointer(
            new com::SocketCommunicationFactory(port, false, network, dir));
        com = comFactory->newCommunication();
#     endif // PRECICE_NO_SOCKETS
    }
    else if (tag.getName() == VALUE_MPI){
      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
#     ifdef PRECICE_NO_MPI
        std::ostringstream error;
        error << "Communication type \"" << VALUE_MPI << "\" can only be used "
              << "when preCICE is compiled with argument \"mpi=on\"";
        throw error.str();
#     else
      comFactory = com::CommunicationFactory::SharedPointer(
          new com::MPIPortsCommunicationFactory(dir));
      com = comFactory->newCommunication();
#     endif
    }
    else if (tag.getName() == VALUE_MPI_SINGLE){
#     ifdef PRECICE_NO_MPI
        std::ostringstream error;
        error << "Communication type \"" << VALUE_MPI_SINGLE << "\" can only be used "
              << "when preCICE is compiled with argument \"mpi=on\"";
        throw error.str();
#     else
        com = com::Communication::SharedPointer(new com::MPIDirectCommunication());
#     endif
    }
    else if (tag.getName() == VALUE_FILES){
      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
      com = com::Communication::SharedPointer(new com::FileCommunication(false, dir));
    }

    assertion(com.get() != NULL);


    DistributedComFactory::SharedPointer distrFactory;
    if(tag.getName() == VALUE_MPI_SINGLE || tag.getName() == VALUE_FILES || distrType == VALUE_GATHER_SCATTER){
      assertion(distrType == VALUE_GATHER_SCATTER);
      distrFactory = DistributedComFactory::SharedPointer(new GatherScatterComFactory(com));
    }
    else if(distrType == VALUE_POINT_TO_POINT){
      assertion(tag.getName() == VALUE_MPI || tag.getName() == VALUE_SOCKETS);
      distrFactory = DistributedComFactory::SharedPointer(new PointToPointComFactory(comFactory));
    }
    assertion(distrFactory.get() != NULL);

    m2n::M2N::SharedPointer m2n = m2n::M2N::SharedPointer(new m2n::M2N(com, distrFactory));
    _m2ns.push_back(boost::make_tuple(m2n, from, to));
  }
}

void M2NConfiguration:: checkDuplicates
(
  const std::string& from,
  const std::string& to )
{
  using boost::get;
  bool alreadyAdded = false;
  for (M2NTuple& tuple : _m2ns){
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
