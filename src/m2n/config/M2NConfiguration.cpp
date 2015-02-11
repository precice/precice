// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "M2NConfiguration.hpp"
#include "m2n/M2N.hpp"
#include "m2n/DistributedComFactory.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "com/SocketCommunicationFactory.hpp"
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
    doc = "Port number to be used by for socket communiation. TODO";
    attrPort.setDocumentation(doc);
    attrPort.setDefaultValue(51235);
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


  XMLAttribute<std::string> attrDistrType ( ATTR_DISTRIBUTION_TYPE);
  doc = "Distribution manner of the M2N communication .";
  attrDistrType.setDocumentation(doc);
  ValidatorEquals<std::string> validDistrGatherScatter ( VALUE_GATHER_SCATTER );
  ValidatorEquals<std::string> validDistrP2P ( VALUE_POINT_TO_POINT);
  attrDistrType.setValidator ( validDistrGatherScatter || validDistrP2P );
  attrDistrType.setDefaultValue(VALUE_POINT_TO_POINT);
  XMLAttribute<std::string> attrFrom ( ATTR_FROM );
  doc = "First participant name involved in communication.";
  attrFrom.setDocumentation(doc);
  XMLAttribute<std::string> attrTo(ATTR_TO);
  doc = "Second participant name involved in communication.";
  attrTo.setDocumentation(doc);

  foreach (XMLTag& tag, tags){
    tag.addAttribute(attrFrom);
    tag.addAttribute(attrTo);
    tag.addAttribute(attrDistrType);
    parent.addSubtag(tag);
  }
}

m2n::PtrM2N M2NConfiguration:: getM2N
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

    com::PtrCommunicationFactory comFactory;
    if (tag.getName() == VALUE_SOCKETS){
#     ifdef PRECICE_NO_SOCKETS
        std::ostringstream error;
        error << "M2N type \"" << VALUE_SOCKETS << "\" can only be used "
              << "when preCICE is compiled with argument \"sockets=on\"";
        throw error.str();
#     else
        std::string network = tag.getStringAttributeValue(ATTR_NETWORK);
        int port = tag.getIntAttributeValue(ATTR_PORT);
        std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
        comFactory = com::PtrCommunicationFactory(new com::SocketCommunicationFactory(network, port, dir));
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
        //TODO
        //comFactory = com::PtrCommunicationFactory(new com::MPIPortsCommunicationFactory(dir));
#     endif
    }
    assertion(comFactory.get() != NULL);

    com::PtrCommunication com = comFactory->newCommunication();

    PtrDistributedComFactory distrFactory;
    if(distrType == VALUE_GATHER_SCATTER){
      distrFactory = PtrDistributedComFactory(new GatherScatterComFactory(com));
    }
    else if(distrType == VALUE_POINT_TO_POINT){
      //TODO
//      distrFactory = PtrDistributedCommunicationFactory(new PointToPointComFactory(comFactory));
    }
    assertion(distrFactory.get() != NULL);

    m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(com, distrFactory));
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
