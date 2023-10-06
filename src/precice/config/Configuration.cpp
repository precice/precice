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
#include "precice/config/SharedPointer.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/ParticipantState.hpp"
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

  auto attrMinTimeStepSize = xml::makeXMLAttribute("min-time-step-size", math::NUMERICAL_ZERO_DIFFERENCE)
                                 .setDocumentation("The smallest maximal time step that preCICE will return."
                                                   "This means that preCICE will round up the time window end by min-time-step-size, which will  "
                                                   " introduce an extra coupling error. For more details see https://github.com/precice/precice/issues/1832.");
  _tag.addAttribute(attrMinTimeStepSize);

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
    _participantConfiguration->setExperimental(_experimental);

    _minTimeStepSize = tag.getDoubleAttributeValue("min-time-step-size");
    PRECICE_CHECK(_minTimeStepSize >= math::NUMERICAL_ZERO_DIFFERENCE, "The minimal time step has to be larger or equal to {}. Please adjust the tag min-time-step-size in the config file.", math::NUMERICAL_ZERO_DIFFERENCE);
    _couplingSchemeConfiguration->setMinTimeStepSize(_minTimeStepSize);

  } else {
    PRECICE_UNREACHABLE("Received callback from unknown tag '{}'.", tag.getName());
  }
}

void Configuration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getName());
  if (tag.getName() == "precice-configuration") {
    //test if both participants do have the exchange meshes
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
  }
}

const PtrParticipantConfiguration &
Configuration::getParticipantConfiguration() const
{
  return _participantConfiguration;
}

} // namespace precice::config
