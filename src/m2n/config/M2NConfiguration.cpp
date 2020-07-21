#include "M2NConfiguration.hpp"
#include <list>
#include <ostream>
#include <stdexcept>
#include "com/CommunicationFactory.hpp"
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/MPISinglePortsCommunicationFactory.hpp"
#include "com/SharedPointer.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/DistributedComFactory.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/M2N.hpp"
#include "m2n/PointToPointComFactory.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"
#include "utils/networking.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

#ifndef PRECICE_NO_MPI
#include <mpi.h>
#endif

namespace precice {
namespace m2n {
M2NConfiguration::M2NConfiguration(xml::XMLTag &parent)
{
  using namespace xml;
  std::string        doc;
  std::list<XMLTag>  tags;
  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  {
    XMLTag tag(*this, "sockets", occ, TAG);
    doc = "Communication via Sockets.";
    tag.setDocumentation(doc);

    auto attrPort = makeXMLAttribute("port", 0)
                        .setDocumentation(
                            "Port number (16-bit unsigned integer) to be used for socket "
                            "communication. The default is \"0\", what means that the OS will "
                            "dynamically search for a free port (if at least one exists) and "
                            "bind it automatically.");
    tag.addAttribute(attrPort);

    auto attrNetwork = makeXMLAttribute("network", utils::networking::loopbackInterfaceName())
                           .setDocumentation(
                               "Interface name to be used for socket communiation. "
                               "Default is the cannonical name of the loopback interface of your platform. "
                               "Might be different on supercomputing systems, e.g. \"ib0\" "
                               "for the InfiniBand on SuperMUC. ");
    tag.addAttribute(attrNetwork);

    auto attrExchangeDirectory = makeXMLAttribute(ATTR_EXCHANGE_DIRECTORY, "")
                                     .setDocumentation(
                                         "Directory where connection information is exchanged. By default, the "
                                         "directory of startup is chosen, and both solvers have to be started "
                                         "in the same directory.");
    tag.addAttribute(attrExchangeDirectory);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, "mpi", occ, TAG);
    doc = "Communication via MPI with startup in separated communication spaces, using multiple communicators.";
    tag.setDocumentation(doc);

    auto attrExchangeDirectory = makeXMLAttribute(ATTR_EXCHANGE_DIRECTORY, "")
                                     .setDocumentation(
                                         "Directory where connection information is exchanged. By default, the "
                                         "directory of startup is chosen, and both solvers have to be started "
                                         "in the same directory.");
    tag.addAttribute(attrExchangeDirectory);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, "mpi-singleports", occ, TAG);
    doc = "Communication via MPI with startup in separated communication spaces, using a single communicator";
    tag.setDocumentation(doc);

    auto attrExchangeDirectory = makeXMLAttribute(ATTR_EXCHANGE_DIRECTORY, "")
                                     .setDocumentation(
                                         "Directory where connection information is exchanged. By default, the "
                                         "directory of startup is chosen, and both solvers have to be started "
                                         "in the same directory.");
    tag.addAttribute(attrExchangeDirectory);
    tags.push_back(tag);
  }

  XMLAttribute<bool> attrEnforce(ATTR_ENFORCE_GATHER_SCATTER, false);
  attrEnforce.setDocumentation("Enforce the distributed communication to a gather-scatter scheme. "
                               "Only recommended for trouble shooting.");

  XMLAttribute<bool> attrTwoLevel(ATTR_USE_TWO_LEVEL_INIT, false);
  attrTwoLevel.setDocumentation("Use a two-level initialization scheme. "
                                "Recommended for large parallel runs (>5000 MPI ranks).");

  auto attrFrom = XMLAttribute<std::string>("from")
                      .setDocumentation(
                          "First participant name involved in communication. For performance reasons, we recommend to use "
                          "the participant with less ranks at the coupling interface as \"from\" in the m2n communication.");
  auto attrTo = XMLAttribute<std::string>("to")
                    .setDocumentation("Second participant name involved in communication.");

  for (XMLTag &tag : tags) {
    tag.addAttribute(attrFrom);
    tag.addAttribute(attrTo);
    tag.addAttribute(attrEnforce);
    tag.addAttribute(attrTwoLevel);
    parent.addSubtag(tag);
  }
}

m2n::PtrM2N M2NConfiguration::getM2N(const std::string &from, const std::string &to)
{
  using std::get;
  for (M2NTuple &tuple : _m2ns) {
    if ((get<1>(tuple) == from) && (get<2>(tuple) == to)) {
      return get<0>(tuple);
    } else if ((get<2>(tuple) == from) && (get<1>(tuple) == to)) {
      return get<0>(tuple);
    }
  }
  PRECICE_ERROR("There is no m2n communication configured between participants \"" + from + "\" and \"" + to + "\". Please add an appropriate \"<m2n />\" tag.");
}

void M2NConfiguration::xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag)
{
  if (tag.getNamespace() == TAG) {
    std::string from = tag.getStringAttributeValue("from");
    std::string to   = tag.getStringAttributeValue("to");
    checkDuplicates(from, to);
    bool enforceGatherScatter = tag.getBooleanAttributeValue(ATTR_ENFORCE_GATHER_SCATTER);
    bool useTwoLevelInit      = tag.getBooleanAttributeValue(ATTR_USE_TWO_LEVEL_INIT);

    if (enforceGatherScatter && useTwoLevelInit) {
      throw std::runtime_error{std::string{"A gather-scatter m2n communication cannot use two-level initialization. Please switch either "} + "\"" + ATTR_ENFORCE_GATHER_SCATTER + "\" or \"" + ATTR_USE_TWO_LEVEL_INIT + "\" off."};
    }
    if (context.size == 1 && useTwoLevelInit) {
      throw std::runtime_error{"To use two-level initialization, both participants need to run in parallel. If you want to run in serial please switch two-level intialization off."};
    }

    com::PtrCommunicationFactory comFactory;
    com::PtrCommunication        com;
    if (tag.getName() == "sockets") {
      std::string network = tag.getStringAttributeValue("network");
      int         port    = tag.getIntAttributeValue("port");

      PRECICE_CHECK(not utils::isTruncated<unsigned short>(port),
                    "The value given for the \"port\" attribute is not a 16-bit unsigned integer: " << port);

      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
      comFactory      = std::make_shared<com::SocketCommunicationFactory>(port, false, network, dir);
      com             = comFactory->newCommunication();
    } else if (tag.getName() == "mpi") {
      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
#ifdef PRECICE_NO_MPI
      PRECICE_ERROR("Communication type \"mpi\" can only be used if preCICE was compiled with MPI support enabled. "
                    "Either switch to a \"sockets\" communication or recompile preCICE with \"PRECICE_MPICommunication=ON\".");
#else
#ifdef OMPI_MAJOR_VERSION
      PRECICE_WARN("preCICE was compiled with OpenMPI and configured to use <m2n:mpi />, which can cause issues in connection build-up. Consider switching to sockets if you encounter problems.");
#endif
      comFactory = std::make_shared<com::MPIPortsCommunicationFactory>(dir);
      com        = comFactory->newCommunication();
#endif
    } else if (tag.getName() == "mpi-singleports") {
      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
#ifdef PRECICE_NO_MPI
      PRECICE_ERROR("Communication type \"mpi-singleports\" can only be used if preCICE was compiled with MPI support enabled. "
                    "Either switch to a \"sockets\" communication or recompile preCICE with \"PRECICE_MPICommunication=ON\".");
#else
#ifdef OMPI_MAJOR_VERSION
      PRECICE_WARN("preCICE was compiled with OpenMPI and configured to use <m2n:mpi-singleports />, which can cause issues in connection build-up. Consider switching to sockets if you encounter problems.");
#endif
      comFactory = std::make_shared<com::MPISinglePortsCommunicationFactory>(dir);
      com        = comFactory->newCommunication();
#endif
    }

    PRECICE_ASSERT(com.get() != nullptr);

    DistributedComFactory::SharedPointer distrFactory;
    if (enforceGatherScatter) {
      distrFactory = std::make_shared<GatherScatterComFactory>(com);
    } else {
      distrFactory = std::make_shared<PointToPointComFactory>(comFactory);
    }
    PRECICE_ASSERT(distrFactory.get() != nullptr);

    auto m2n = std::make_shared<m2n::M2N>(com, distrFactory, false, useTwoLevelInit);
    _m2ns.push_back(std::make_tuple(m2n, from, to));
  }
}

void M2NConfiguration::checkDuplicates(
    const std::string &from,
    const std::string &to)
{
  using std::get;
  bool alreadyAdded = false;
  for (M2NTuple &tuple : _m2ns) {
    alreadyAdded |= (get<1>(tuple) == from) && (get<2>(tuple) == to);
    alreadyAdded |= (get<2>(tuple) == from) && (get<1>(tuple) == to);
  }
  PRECICE_CHECK(!alreadyAdded, "Multiple m2n communications between participant \"" + from + "\" and \"" + to + "\" are not allowed. Please remove redundant <m2n /> tags between them.");
}

} // namespace m2n
} // namespace precice
