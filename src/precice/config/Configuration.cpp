#include "Configuration.hpp"
#include <map>
#include <memory>
#include <ostream>
#include <vector>
#include "ParticipantConfiguration.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "math/differences.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "partition/ProvidedPartition.hpp"
#include "partition/ReceivedPartition.hpp"
#include "precice/config/SharedPointer.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/ParticipantState.hpp"
#include "precice/impl/ProvidedMeshContext.hpp"
#include "precice/impl/ReceivedMeshContext.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice::config {

Configuration::Configuration()
    : _tag(*this, "precice-configuration", xml::XMLTag::OCCUR_ONCE),
      _logConfig(_tag), // This must be the first configuration to be constructed
      _profilingConfig(_tag)
{
  _tag.setDocumentation("Main tag containing preCICE configuration.");
  _tag.addNamespace("data");
  _tag.addNamespace("communication");
  _tag.addNamespace("mapping");
  _tag.addNamespace("export");
  _tag.addNamespace("action");
  _tag.addNamespace("coupling-scheme");
  _tag.addNamespace("acceleration");

  auto attrExperimental = xml::makeXMLAttribute("experimental", false)
                              .setDocumentation("Enable experimental features.");
  _tag.addAttribute(attrExperimental);

  auto attrRemeshing = xml::makeXMLAttribute("allow-remeshing", false)
                           .setDocumentation("Enable experimental remeshing feature, requires experimental to be true.");
  _tag.addAttribute(attrRemeshing);

  auto attrWaitInFinalize = xml::makeXMLAttribute("wait-in-finalize", false)
                                .setDocumentation("Connected participants wait for each other in finalize, which can be helpful in SLURM sessions.");
  _tag.addAttribute(attrWaitInFinalize);
  _dataConfiguration = std::make_shared<mesh::DataConfiguration>(
      _tag);
  _meshConfiguration = std::make_shared<mesh::MeshConfiguration>(
      _tag, _dataConfiguration);
  _m2nConfiguration = std::make_shared<m2n::M2NConfiguration>(
      _tag);
  _participantConfiguration = std::make_shared<ParticipantConfiguration>(
      _tag, _meshConfiguration);
  _couplingSchemeConfiguration = std::make_shared<cplscheme::CouplingSchemeConfiguration>(
      _tag, _meshConfiguration, _m2nConfiguration, _participantConfiguration);
}

xml::XMLTag &Configuration::getXMLTag()
{
  return _tag;
}

void Configuration::xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag)
{
  PRECICE_TRACE(tag.getName());
  if (tag.getName() == "precice-configuration") {
    _experimental = tag.getBooleanAttributeValue("experimental");
    _remeshing    = tag.getBooleanAttributeValue("allow-remeshing");

    PRECICE_CHECK(!_remeshing || _experimental, "Remeshing is considered an experimental feature. Please enable <precice-configuration experimental=\"1\" >.");
    _participantConfiguration->setExperimental(_experimental);
    _participantConfiguration->setRemeshing(_remeshing);
    _couplingSchemeConfiguration->setRemeshing(_remeshing);
    _waitInFinalize = tag.getBooleanAttributeValue("wait-in-finalize");
  } else {
    PRECICE_UNREACHABLE("Received callback from unknown tag '{}'.", tag.getName());
  }
}

void Configuration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag                     &tag)
{
  PRECICE_TRACE(tag.getName());
  PRECICE_ASSERT(tag.getName() == "precice-configuration");

  // test if both participants do have the exchange meshes
  typedef std::map<std::string, std::vector<std::string>>::value_type neededMeshPair;
  for (const neededMeshPair &neededMeshes : _meshConfiguration->getNeededMeshes()) {
    bool participantFound = false;
    for (const impl::PtrParticipant &participant : _participantConfiguration->getParticipants()) {
      if (participant->getName() == neededMeshes.first) {
        for (const std::string &neededMesh : neededMeshes.second) {
          PRECICE_CHECK(participant->isMeshUsed(neededMesh),
                        "Participant \"{}\" needs to use the mesh \"{}\" to be able to use it in the coupling scheme. "
                        "Please either add a provide-mesh or a receive-mesh tag in this participant's configuration, or use a different mesh in the coupling scheme.",
                        neededMeshes.first, neededMesh);
        }
        participantFound = true;
        break;
      }
    }
    PRECICE_ASSERT(participantFound);
  }

  // test if all M2Ns use participants that exist
  for (const auto &m2n : _m2nConfiguration->m2ns()) {
    PRECICE_CHECK(_participantConfiguration->hasParticipant(m2n.acceptor),
                  "The acceptor in <m2n:... acceptor=\"{}\" connector=\"{}\" /> is an unknown. {}",
                  m2n.acceptor, m2n.connector, _participantConfiguration->hintFor(m2n.acceptor));

    PRECICE_CHECK(_participantConfiguration->hasParticipant(m2n.connector),
                  "The connector in <m2n:... acceptor=\"{}\" connector=\"{}\" /> is an unknown. {}",
                  m2n.acceptor, m2n.connector, _participantConfiguration->hintFor(m2n.connector));
  }
}

const PtrParticipantConfiguration &
Configuration::getParticipantConfiguration() const
{
  return _participantConfiguration;
}

std::map<std::string, m2n::BoundM2N> Configuration::getBoundM2NsFor(std::string_view participantName) const
{
  std::map<std::string, m2n::BoundM2N> result;

  for (const auto &m2nConf : _m2nConfiguration->m2ns()) {
    if (m2nConf.acceptor != participantName && m2nConf.connector != participantName) {
      continue;
    }

    std::string comPartner("");
    bool        isRequesting;
    if (m2nConf.acceptor == participantName) {
      comPartner   = m2nConf.connector;
      isRequesting = true;
    } else {
      comPartner   = m2nConf.acceptor;
      isRequesting = false;
    }

    PRECICE_ASSERT(!comPartner.empty());
    for (const impl::PtrParticipant &participant : _participantConfiguration->getParticipants()) {
      if (participant->getName() == comPartner) {
        PRECICE_ASSERT(not utils::contained(comPartner, result), comPartner);
        PRECICE_ASSERT(m2nConf.m2n);

        result[comPartner] = [&] {
          m2n::BoundM2N bound;
          bound.m2n          = m2nConf.m2n;
          bound.localName    = participantName;
          bound.remoteName   = comPartner;
          bound.isRequesting = isRequesting;
          return bound;
        }();
      }
    }
  }
  return result;
}

void Configuration::configurePartitionsFor(std::string_view participantName)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_participantConfiguration->hasParticipant(participantName));

  auto participant = _participantConfiguration->getParticipant(participantName);

  // Iterate over provided mesh contexts
  for (auto &context : participant->providedMeshContexts()) {
    // Accessor provides mesh
    context.partition = std::make_shared<partition::ProvidedPartition>(context.mesh);

    for (auto &receiver : _participantConfiguration->getParticipants()) {
      for (auto &receivedMeshContext : receiver->receivedMeshContexts()) {
        if (receivedMeshContext.receiveMeshFrom == participantName && receivedMeshContext.mesh->getName() == context.mesh->getName()) {
          // meshRequirement has to be copied from "from" to provide", since
          // mapping are only defined at "provide"
          if (receivedMeshContext.meshRequirement > context.meshRequirement) {
            context.meshRequirement = receivedMeshContext.meshRequirement;
          }

          m2n::PtrM2N m2n = _m2nConfiguration->getM2N(receiver->getName(), std::string(participantName));
          m2n->createDistributedCommunication(context.mesh);
          context.partition->addM2N(m2n);
        }
      }
    }
  }

  // Iterate over received mesh contexts
  for (auto &receivedContext : participant->receivedMeshContexts()) {
    // Accessor receives mesh
    std::string receiver(participantName);
    std::string provider(receivedContext.receiveMeshFrom);

    PRECICE_DEBUG("Receiving mesh from {}", provider);

    receivedContext.partition = std::make_shared<partition::ReceivedPartition>(receivedContext.mesh, receivedContext.geoFilter, receivedContext.safetyFactor, receivedContext.allowDirectAccess);

    m2n::PtrM2N m2n = _m2nConfiguration->getM2N(receiver, provider);
    m2n->createDistributedCommunication(receivedContext.mesh);
    receivedContext.partition->addM2N(m2n);
    for (const precice::impl::MappingContext &mappingContext : receivedContext.fromMappingContexts) {
      receivedContext.partition->addFromMapping(mappingContext.mapping);
    }
    for (const precice::impl::MappingContext &mappingContext : receivedContext.toMappingContexts) {
      receivedContext.partition->addToMapping(mappingContext.mapping);
    }
  }
}

} // namespace precice::config
