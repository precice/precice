#include "LogConfiguration.hpp"
#include "xml/ValidatorEquals.hpp"
#include "xml/ValidatorOr.hpp"

namespace precice {
namespace config {


LogConfiguration::LogConfiguration
(
  xml::XMLTag& parent)
{
  // We do default initialization here, so logging will be initialized
  // as soon as possible and also if there is no <log> tag.
  precice::logging::setupLogging();
  
  using namespace xml;
  XMLTag tagLog(*this, "log", XMLTag::OCCUR_NOT_OR_ONCE);
  tagLog.setDocumentation("Configures logging");

  XMLAttribute<bool> attrLogEnabled("enabled");
  attrLogEnabled.setDocumentation("Enables logging");
  attrLogEnabled.setDefaultValue(true);
  tagLog.addAttribute(attrLogEnabled);
  
  XMLTag tagSink(*this, "sink", XMLTag::OCCUR_ARBITRARY);
  XMLAttribute<std::string> attrType("type");
  attrType.setDocumentation("Type of sink");
  attrType.setValidator ( ValidatorEquals<std::string>("stream") || ValidatorEquals<std::string>("file") );
  attrType.setDefaultValue(precice::logging::BackendConfiguration::default_type);
  tagSink.addAttribute(attrType);

  XMLAttribute<std::string> attrOutput("output");
  attrOutput.setDocumentation("Output. If type=stream it can be stdout or stderr. Otherwise it is a filename");
  attrOutput.setDefaultValue(precice::logging::BackendConfiguration::default_output);
  tagSink.addAttribute(attrOutput);

  XMLAttribute<std::string> attrFormat("format");
  attrFormat.setDocumentation("Boost Log Format String");
  attrFormat.setDefaultValue(precice::logging::BackendConfiguration::default_formatter);
  tagSink.addAttribute(attrFormat);

  XMLAttribute<std::string> attrFilter("filter");
  attrFilter.setDocumentation("Boost Log Filter String");
  attrFilter.setDefaultValue(precice::logging::BackendConfiguration::default_filter);
  tagSink.addAttribute(attrFilter);
  
  XMLAttribute<bool> attrEnabled("enabled");
  attrEnabled.setDocumentation("Enables the sink");
  attrEnabled.setDefaultValue(true);
  tagSink.addAttribute(attrEnabled);
  
  tagLog.addSubtag(tagSink);
  parent.addSubtag(tagLog);
}

void LogConfiguration::xmlTagCallback
(
  xml::XMLTag& tag )
{
  TRACE(tag.getFullName());
  
  if (tag.getName() == "sink" and tag.getBooleanAttributeValue("enabled")) {
    precice::logging::BackendConfiguration config;
    config.setOption("type", tag.getStringAttributeValue("type"));
    config.setOption("output", tag.getStringAttributeValue("output"));
    config.setOption("filter", tag.getStringAttributeValue("filter"));
    config.setOption("format", tag.getStringAttributeValue("format"));
    config.setOption("enabled", "true"); // Not needed, but correct.
    _logconfig.push_back(config);
  }
}

void LogConfiguration::xmlEndTagCallback
(
  xml::XMLTag& tag )
{
  TRACE(tag.getFullName());
  if (tag.getName() == "log")
    precice::logging::setupLogging(_logconfig, tag.getBooleanAttributeValue("enabled"));
}

}} // namespace precice, config
