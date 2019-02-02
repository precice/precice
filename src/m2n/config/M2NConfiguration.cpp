#include "M2NConfiguration.hpp"
#include <list>
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/MPISinglePortsCommunicationFactory.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "m2n/DistributedComFactory.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/M2N.hpp"
#include "m2n/PointToPointComFactory.hpp"
#include "xml/ValidatorEquals.hpp"
#include "xml/ValidatorOr.hpp"
#include "xml/XMLAttribute.hpp"
#include "utils/Helpers.hpp"

namespace precice
{
namespace m2n
{
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

    XMLAttribute<int> attrPort("port");
    doc = "Port number (16-bit unsigned integer) to be used for socket ";
    doc += "communication. The default is \"0\", what means that the OS will ";
    doc += "dynamically search for a free port (if at least one exists) and ";
    doc += "bind it automatically.";
    attrPort.setDocumentation(doc);
    attrPort.setDefaultValue(0);
    tag.addAttribute(attrPort);

    XMLAttribute<std::string> attrNetwork("network");
    doc = "Interface name to be used for socket communiation. ";
    doc += "Default is \"lo\", i.e., the local host loopback. ";
    doc += "Might be different on supercomputing systems, e.g. \"ib0\" ";
    doc += "for the InfiniBand on SuperMUC. ";
    doc += "For macOS use \"lo0\". ";
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
    XMLTag tag(*this, "mpi", occ, TAG);
    doc = "Communication via MPI with startup in separated communication spaces, using multiple communicators.";
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
    XMLTag tag(*this, "mpi-singleports", occ, TAG);
    doc = "Communication via MPI with startup in separated communication spaces, using a single communicator";
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
    XMLTag tag(*this, "mpi-single", occ, TAG);
    doc = "Communication via MPI with startup in common communication space.";
    tag.setDocumentation(doc);
    tags.push_back(tag);
  }

  XMLAttribute<std::string> attrDistrTypeBoth(ATTR_DISTRIBUTION_TYPE);
  doc = "Distribution manner of the M2N communication. ";
  doc += "\"" + VALUE_POINT_TO_POINT + "\" uses a pure point to point communication and is recommended. ";
  doc += "\"" + VALUE_GATHER_SCATTER + "\" should only be used if at least one serial participant is used ";
  doc += "or for troubleshooting.";
  attrDistrTypeBoth.setDocumentation(doc);
  auto validDistrGatherScatter = makeValidatorEquals(VALUE_GATHER_SCATTER);
  auto validDistrP2P = makeValidatorEquals(VALUE_POINT_TO_POINT);
  attrDistrTypeBoth.setValidator(validDistrGatherScatter || validDistrP2P);
  attrDistrTypeBoth.setDefaultValue(VALUE_POINT_TO_POINT);

  XMLAttribute<std::string> attrDistrTypeOnly(ATTR_DISTRIBUTION_TYPE);
  doc = "Distribution manner of the M2N communication .";
  doc += "\"" + VALUE_POINT_TO_POINT + "\" uses a pure point to point communication and is recommended. ";
  doc += "\"" + VALUE_GATHER_SCATTER + "\" should only be used if at least one serial participant is used ";
  doc += "or for troubleshooting.";
  attrDistrTypeOnly.setDocumentation(doc);
  attrDistrTypeOnly.setValidator(validDistrGatherScatter);
  attrDistrTypeOnly.setDefaultValue(VALUE_GATHER_SCATTER);

  XMLAttribute<std::string> attrFrom("from");
  doc = "First participant name involved in communication. For performance reasons, we recommend to use ";
  doc += "the participant with less ranks at the coupling interface as \"from\" in the m2n communication.";
  attrFrom.setDocumentation(doc);
  XMLAttribute<std::string> attrTo("to");
  doc = "Second participant name involved in communication.";
  attrTo.setDocumentation(doc);

  for (XMLTag &tag : tags) {
    tag.addAttribute(attrFrom);
    tag.addAttribute(attrTo);
    if (tag.getName() == "mpi" || tag.getName() == "mpi-singleports" || tag.getName() == "sockets") {
      tag.addAttribute(attrDistrTypeBoth);
    } else {
      tag.addAttribute(attrDistrTypeOnly);
    }
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
  std::ostringstream error;
  error << "No m2n communication configured between \"" << from << "\" and \"" << to << "\"!";
  throw error.str();
}

void M2NConfiguration::xmlTagCallback(xml::XMLTag &tag)
{
  if (tag.getNamespace() == TAG) {
    std::string from = tag.getStringAttributeValue("from");
    std::string to   = tag.getStringAttributeValue("to");
    checkDuplicates(from, to);
    std::string distrType = tag.getStringAttributeValue(ATTR_DISTRIBUTION_TYPE);

    com::PtrCommunicationFactory comFactory;
    com::PtrCommunication        com;
    if (tag.getName() == "sockets") {
      std::string network = tag.getStringAttributeValue("network");
      int         port    = tag.getIntAttributeValue("port");

      CHECK(not utils::isTruncated<unsigned short>(port),
            "The value given for the \"port\" attribute is not a 16-bit unsigned integer: " << port);

      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
      comFactory      = std::make_shared<com::SocketCommunicationFactory>(port, false, network, dir);
      com             = comFactory->newCommunication();
    } else if (tag.getName() == "mpi") {
      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
#ifdef PRECICE_NO_MPI
      std::ostringstream error;
      error << "Communication type \"mpi\" can only be used"
            << "when preCICE is compiled with argument \"mpi=on\"";
      throw error.str();
#else
      comFactory = std::make_shared<com::MPIPortsCommunicationFactory>(dir);
      com        = comFactory->newCommunication();
#endif
    } else if (tag.getName() == "mpi-singleports") {
      std::string dir = tag.getStringAttributeValue(ATTR_EXCHANGE_DIRECTORY);
#ifdef PRECICE_NO_MPI
      std::ostringstream error;
      error << "Communication type \"mpi-singleports\" can only be used "
            << "when preCICE is compiled with argument \"mpi=on\"";
      throw error.str();
#else
      comFactory = std::make_shared<com::MPISinglePortsCommunicationFactory>(dir);
      com        = comFactory->newCommunication();
#endif
    } else if (tag.getName() == "mpi-single") {
#ifdef PRECICE_NO_MPI
      std::ostringstream error;
      error << "Communication type \"" << "mpi-single" << "\" can only be used "
            << "when preCICE is compiled with argument \"mpi=on\"";
      throw error.str();
#else
      com        = std::make_shared<com::MPIDirectCommunication>();
#endif
    }

    assertion(com.get() != nullptr);

    DistributedComFactory::SharedPointer distrFactory;
    if (tag.getName() == "mpi-single" || distrType == VALUE_GATHER_SCATTER) {
      assertion(distrType == VALUE_GATHER_SCATTER);
      distrFactory = std::make_shared<GatherScatterComFactory>(com);
    } else if (distrType == VALUE_POINT_TO_POINT) {
      assertion(tag.getName() == "mpi" or tag.getName() == "mpi-singleports" or tag.getName() == "sockets");
      distrFactory = std::make_shared<PointToPointComFactory>(comFactory);
    }
    assertion(distrFactory.get() != nullptr);

    auto m2n = std::make_shared<m2n::M2N>(com, distrFactory);
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
  if (alreadyAdded) {
    std::ostringstream error;
    error << "Multiple communication defined between participant \"" << from
          << "\" and \"" << to << "\"";
    throw error.str();
  }
}

} // namespace m2n
} // namespace precice
