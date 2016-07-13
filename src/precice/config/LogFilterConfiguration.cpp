#include "LogFilterConfiguration.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"

namespace precice {
namespace config {

tarch::logging::Log LogFilterConfiguration::
  _log("precice::config::LogFilterConfiguration");

LogFilterConfiguration:: LogFilterConfiguration
(
  utils::XMLTag& parent )
:
  TAG("log-filter"),
  ATTR_TARGET("target"),
  ATTR_COMPONENT("component"),
  ATTR_SWITCH("switch"),
  _filterList()
{
  using namespace utils;
  std::string doc;
  XMLTag tag(*this, TAG, XMLTag::OCCUR_ONCE_OR_MORE);
  doc = "Filters debug and/or info messages written by preCICE. ";
  doc += "Debug messages are only written when preCICE is compiled in debug mode. ";
  doc += "There needs to be at least one filter entry valid for any output ";
  doc += "message. This can be achieved by using the empty string for ";
  doc += "attributes \"target\" and \"component\".";
  tag.setDocumentation(doc);

  XMLAttribute<std::string> attrTarget(ATTR_TARGET);
  doc = "Specifies the message type/s to be filtered. ";
  doc += "An empty string means debug and info.";
  attrTarget.setDocumentation(doc);
  ValidatorEquals<std::string> validDebug("debug");
  ValidatorEquals<std::string> validInfo("info");
  ValidatorEquals<std::string> validAll("");
  attrTarget.setValidator(validDebug || validInfo || validAll);
  tag.addAttribute(attrTarget);

  XMLAttribute<std::string> attrComponent(ATTR_COMPONENT);
  doc = "Specifies the component to be filtered, which is visible when log-trace";
  doc += " in the <log-output> tag is set to true. It is possible to specify ";
  doc += "only a part of the full component name (e.g. just the namespace) to ";
  doc += "address all components starting with that part. The empty string means";
  doc += "all components.";
  attrComponent.setDocumentation(doc);
  tag.addAttribute(attrComponent);

  XMLAttribute<bool> attrSwitch(ATTR_SWITCH);
  doc = "If switch is set to \"on\", the messages are enabled, if set to \"off\"";
  doc += " messages are disabled. By default, all messages are disabled.";
  attrSwitch.setDocumentation(doc);
  tag.addAttribute(attrSwitch);

  parent.addSubtag(tag);
}

tarch::logging::CommandLineLogger::FilterList LogFilterConfiguration:: getFilterList() const
{
  return _filterList;
}

void LogFilterConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlTagCallback()", tag.getFullName());
  if (tag.getName() == TAG){
    tarch::logging::CommandLineLogger::FilterListEntry newEntry;
    newEntry._targetName = tag.getStringAttributeValue(ATTR_TARGET);
    newEntry._isBlackEntry = not tag.getBooleanAttributeValue(ATTR_SWITCH);
    newEntry._rank = -1; // no ranks configured
    newEntry._namespaceName = tag.getStringAttributeValue(ATTR_COMPONENT);
    _filterList.insert(newEntry);
  }
}

void LogFilterConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlEndTagCallback()", tag.getFullName());
}

}} // namespace precice, config
