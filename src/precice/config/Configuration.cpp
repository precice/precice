#include "Configuration.hpp"
#include "logging/LogMacros.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice {

extern bool syncMode;

namespace config {

Configuration::Configuration()
    : _tag(*this, "precice-configuration", xml::XMLTag::OCCUR_ONCE),
      _logConfig(_tag),
      _solverInterfaceConfig(_tag)
{
  _tag.setDocumentation("Main tag containing preCICE configuration.");
  _tag.addNamespace("data");
  _tag.addNamespace("communication");
  _tag.addNamespace("mapping");
  _tag.addNamespace("export");
  _tag.addNamespace("action");
  _tag.addNamespace("coupling-scheme");
  _tag.addNamespace("acceleration");

  auto attrSyncMode = xml::makeXMLAttribute("sync-mode", false)
                          .setDocumentation("sync-mode enabled additional inter- and intra-participant synchronizations");
  _tag.addAttribute(attrSyncMode);
}

xml::XMLTag &Configuration::getXMLTag()
{
  return _tag;
}

void Configuration::xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag)
{
  PRECICE_TRACE(tag.getName());
  if (tag.getName() == "precice-configuration") {
    precice::syncMode = tag.getBooleanAttributeValue("sync-mode");
  }
}

void Configuration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getName());
}

const SolverInterfaceConfiguration &
Configuration::getSolverInterfaceConfiguration() const
{
  return _solverInterfaceConfig;
}

} // namespace config
} // namespace precice
