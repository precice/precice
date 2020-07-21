#include "LogConfiguration.hpp"
#include "logging/LogMacros.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice {
namespace config {

LogConfiguration::LogConfiguration(
    xml::XMLTag &parent)
{
  // We do default initialization here, so logging will be initialized
  // as soon as possible and also if there is no <log> tag.
  precice::logging::setupLogging();

  using namespace xml;
  XMLTag tagLog(*this, "log", XMLTag::OCCUR_NOT_OR_ONCE);
  tagLog.setDocumentation("Configures logging");

  auto attrLogEnabled = makeXMLAttribute("enabled", true)
                            .setDocumentation("Enables logging");
  tagLog.addAttribute(attrLogEnabled);

  XMLTag tagSink(*this, "sink", XMLTag::OCCUR_ARBITRARY);
  auto   attrType = XMLAttribute<std::string>("type")
                      .setDocumentation("Type of sink")
                      .setOptions({"stream", "file"})
                      .setDefaultValue(precice::logging::BackendConfiguration::default_type);
  tagSink.addAttribute(attrType);

  auto attrOutput = XMLAttribute<std::string>("output")
                        .setDocumentation("Output. If type=stream it can be stdout or stderr. Otherwise it is a filename")
                        .setDefaultValue(precice::logging::BackendConfiguration::default_output);
  tagSink.addAttribute(attrOutput);

  auto attrFormat = XMLAttribute<std::string>("format")
                        .setDocumentation("Boost Log Format String")
                        .setDefaultValue(precice::logging::BackendConfiguration::default_formatter);
  tagSink.addAttribute(attrFormat);

  auto attrFilter = XMLAttribute<std::string>("filter")
                        .setDocumentation("Boost Log Filter String")
                        .setDefaultValue(precice::logging::BackendConfiguration::default_filter);
  tagSink.addAttribute(attrFilter);

  auto attrEnabled = makeXMLAttribute("enabled", true)
                         .setDocumentation("Enables the sink");
  tagSink.addAttribute(attrEnabled);

  tagLog.addSubtag(tagSink);
  parent.addSubtag(tagLog);
}

void LogConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getFullName());

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

void LogConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getFullName());
  if (tag.getName() == "log")
    precice::logging::setupLogging(_logconfig, tag.getBooleanAttributeValue("enabled"));
}

} // namespace config
} // namespace precice
