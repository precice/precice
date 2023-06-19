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

namespace precice::m2n {
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
                               "Interface name to be used for socket communication. "
                               "Default is the canonical name of the loopback interface of your platform. "
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
    XMLTag tag(*this, "mpi-multiple-ports", occ, TAG);
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
    XMLTag tag(*this, "mpi", occ, TAG);
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
  {
    /// @TODO Remove in Version 3.0, see https://github.com/precice/precice/issues/1650
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

  auto attrFrom = XMLAttribute<std::string>("acceptor")
                      .setDocumentation(
                          "First participant name involved in communication. For performance reasons, we recommend to use "
                          "the participant with less ranks at the coupling interface as \"acceptor\" in the m2n communication.");
  auto attrTo = XMLAttribute<std::string>("connector")
                    .setDocumentation("Second participant name involved in communication.");

  for (XMLTag &tag : tags) {
    tag.addAttribute(attrFrom);
    tag.addAttribute(attrTo);
    tag.addAttribute(attrEnforce);
    tag.addAttribute(attrTwoLevel);
    parent.addSubtag(tag);
  }
}

m2n::PtrM2N M2NConfiguration::getM2N(const std::string &acceptor, const std::string &connector)
{
  using std::get;
  for (M2NTuple &tuple : _m2ns) {
    if ((get<1>(tuple) == acceptor) && (get<2>(tuple) == connector)) {
      return get<0>(tuple);
    } else if ((get<2>(tuple) == acceptor) && (get<1>(tuple) == connector)) {
      return get<0>(tuple);
    }
  }
  PRECICE_ERROR("There is no m2n communication configured between participants \"" + acceptor + "\" and \"" + connector + "\". Please add an appropriate \"<m2n />\" tag.");
}

bool M2NConfiguration::isM2NConfigured(const std::string &acceptor, const std::string &connector)
{
  return std::any_of(std::begin(_m2ns), std::end(_m2ns),
                     [acceptor, connector](const auto &m2nTuple) {
                       return ((std::get<1>(m2nTuple) == acceptor) && (std::get<2>(m2nTuple) == connector)) || ((std::get<1>(m2nTuple) == connector) && (std::get<2>(m2nTuple) == acceptor));
                     });
}

void M2NConfiguration::xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag)
{
  if (tag.getNamespace() == TAG) {
    std::string acceptor  = tag.getStringAttributeValue("acceptor");
    std::string connector = tag.getStringAttributeValue("connector");
    checkDuplicates(acceptor, connector);
    bool enforceGatherScatter = tag.getBooleanAttributeValue(ATTR_ENFORCE_GATHER_SCATTER);
    bool useTwoLevelInit      = tag.getBooleanAttributeValue(ATTR_USE_TWO_LEVEL_INIT);

    if (enforceGatherScatter && useTwoLevelInit) {
      throw std::runtime_error{std::string{"A gather-scatter m2n communication cannot use two-level initialization. Please switch either "} + "\"" + ATTR_ENFORCE_GATHER_SCATTER + "\" or \"" + ATTR_USE_TWO_LEVEL_INIT + "\" off."};
    }
    if (context.size == 1 && useTwoLevelInit) {
      throw std::runtime_error{"To use two-level initialization, both participants need to run in parallel. If you want to run in serial please switch two-level initialization off."};
    }

    com::PtrCommunicationFactory comFactory;
    com::PtrCommunication        com;
    const std::string            tagName = tag.getName();
    if (tagName == "sockets") {
      std::string network = tag.getStringAttributeValue("network");
      int         port    = tag.getIntAttributeValue("port");

      PRECICE_CHECK(not utils::isTruncated<unsigned short>(port),
                    "The value given for the \"port\" attribute is not a 16-bit unsigned integer: {}", port);

      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
      comFactory      = std::make_shared<com::SocketCommunicationFactory>(port, false, network, dir);
      com             = comFactory->newCommunication();
    } else if (tagName == "mpi-multiple-ports") {
      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
#ifdef PRECICE_NO_MPI
      PRECICE_ERROR("Communication type \"mpi-multiple-ports\" can only be used if preCICE was compiled with MPI support enabled. "
                    "Either switch to a \"sockets\" communication or recompile preCICE with \"PRECICE_MPICommunication=ON\".");
#else
#ifdef OMPI_MAJOR_VERSION
      PRECICE_WARN("preCICE was compiled with OpenMPI and configured to use <m2n:mpi-multiple-ports />, which can cause issues in connection build-up. Consider switching to sockets if you encounter problems. Ignore this warning if participants find each other and the simulation starts.");
#endif
      comFactory = std::make_shared<com::MPIPortsCommunicationFactory>(dir);
      com        = comFactory->newCommunication();
#endif
    } else if (tagName == "mpi" || tagName == "mpi-singleports") {
      if (tagName == "mpi-singleports") {
        PRECICE_WARN("You used <m2n:mpi-singleports />, which is deprecated. Please use <m2n:mpi /> instead.");
      }
      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
#ifdef PRECICE_NO_MPI
      PRECICE_ERROR("Communication type \"{}\" can only be used if preCICE was compiled with MPI support enabled. "
                    "Either switch to a \"sockets\" communication or recompile preCICE with \"PRECICE_MPICommunication=ON\".",
                    tagName);
#else
#ifdef OMPI_MAJOR_VERSION
      PRECICE_WARN("preCICE was compiled with OpenMPI and configured to use <m2n:{} />, which can cause issues in connection build-up. Consider switching to sockets if you encounter problems. Ignore this warning if participants find each other and the simulation starts.", tagName);
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
    _m2ns.emplace_back(m2n, acceptor, connector);
  }
}

void M2NConfiguration::checkDuplicates(
    const std::string &acceptor,
    const std::string &connector)
{
  using std::get;
  bool alreadyAdded = false;
  for (M2NTuple &tuple : _m2ns) {
    alreadyAdded |= (get<1>(tuple) == acceptor) && (get<2>(tuple) == connector);
    alreadyAdded |= (get<2>(tuple) == acceptor) && (get<1>(tuple) == connector);
  }
  PRECICE_CHECK(!alreadyAdded, "Multiple m2n communications between participant \"" + acceptor + "\" and \"" + connector + "\" are not allowed. Please remove redundant <m2n /> tags between them.");
}

} // namespace precice::m2n
