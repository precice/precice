#include "SolverInterfaceConfiguration.hpp"
#include <map>
#include <memory>
#include <ostream>
#include <vector>
#include "ParticipantConfiguration.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "precice/config/SharedPointer.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice::config {

SolverInterfaceConfiguration::SolverInterfaceConfiguration(xml::XMLTag &parent)
{
  using namespace xml;
  XMLTag tag(*this, "solver-interface", XMLTag::OCCUR_ONCE);
  tag.setDocumentation("Configuration of simulation relevant features.");
  auto attrDimensions = makeXMLAttribute("dimensions", 2)
                            .setDocumentation("Determines the spatial dimensionality of the configuration")
                            .setOptions({2, 3});
  tag.addAttribute(attrDimensions);
  auto attrExperimental = makeXMLAttribute("experimental", false)
                              .setDocumentation("Enable experimental features.");
  tag.addAttribute(attrExperimental);

  _dataConfiguration = std::make_shared<mesh::DataConfiguration>(
      tag);
  _meshConfiguration = std::make_shared<mesh::MeshConfiguration>(
      tag, _dataConfiguration);
  _m2nConfiguration = std::make_shared<m2n::M2NConfiguration>(
      tag);
  _participantConfiguration = std::make_shared<ParticipantConfiguration>(
      tag, _meshConfiguration);
  _couplingSchemeConfiguration = std::make_shared<cplscheme::CouplingSchemeConfiguration>(
      tag, _meshConfiguration, _m2nConfiguration, _participantConfiguration);

  parent.addSubtag(tag);
}

void SolverInterfaceConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE();
  if (tag.getName() == "solver-interface") {
    _dimensions = tag.getIntAttributeValue("dimensions");
    _dataConfiguration->setDimensions(_dimensions);
    _meshConfiguration->setDimensions(_dimensions);
    _participantConfiguration->setDimensions(_dimensions);
    _experimental = tag.getBooleanAttributeValue("experimental");
    _couplingSchemeConfiguration->setExperimental(_experimental);
    _participantConfiguration->setExperimental(_experimental);
  } else {
    PRECICE_UNREACHABLE("Received callback from unknown tag '{}'.", tag.getName());
  }
}

void SolverInterfaceConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE();
  if (tag.getName() == "solver-interface") {
    //test if both participants do have the exchange meshes
    typedef std::map<std::string, std::vector<std::string>>::value_type neededMeshPair;
    for (const neededMeshPair &neededMeshes : _meshConfiguration->getNeededMeshes()) {
      bool participantFound = false;
      for (const impl::PtrParticipant &participant : _participantConfiguration->getParticipants()) {
        if (participant->getName() == neededMeshes.first) {
          for (const std::string &neededMesh : neededMeshes.second) {
            PRECICE_CHECK(participant->isMeshUsed(neededMesh),
                          ::precice::ConfigurationError,
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
  }
}

int SolverInterfaceConfiguration::getDimensions() const
{
  return _dimensions;
}

const PtrParticipantConfiguration &
SolverInterfaceConfiguration::getParticipantConfiguration() const
{
  return _participantConfiguration;
}

} // namespace precice::config
