#include "LogConfiguration.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"

#include <iostream>

namespace precice {
namespace config {

/// Logging device.
precice::logging::Logger precice::config::LogConfiguration::_log("logging::config::LogConfiguration");

LogConfiguration::LogConfiguration
(
  utils::XMLTag& parent)
{
  // We do default initialization here, so logging will be initialized
  // as soon as possible and also if there is no <log> tag.
  precice::logging::setupLogging();
  
  using namespace utils;
  XMLTag tagLog(*this, "log", XMLTag::OCCUR_NOT_OR_ONCE);
  tagLog.setDocumentation("Configures logging");

  XMLTag tagSink(*this, "sink", XMLTag::OCCUR_ARBITRARY);
  XMLAttribute<std::string> attrType("type");
  attrType.setDocumentation("Type of sink");
  attrType.setValidator ( ValidatorEquals<std::string>("stream") || ValidatorEquals<std::string>("file") );
  tagSink.addAttribute(attrType);

  XMLAttribute<std::string> attrOutput("output");
  attrOutput.setDocumentation("Output. If type=stream it can be stdout or stderr. Otherwise it is a filename");
  tagSink.addAttribute(attrOutput);

  XMLAttribute<std::string> attrFormat("format");
  attrFormat.setDocumentation("Boost Log Format String");
  attrFormat.setDefaultValue(precice::logging::BackendConfiguration::default_formatter);
  tagSink.addAttribute(attrFormat);

  XMLAttribute<std::string> attrFilter("filter");
  attrFilter.setDocumentation("Boost Log Filter String");
  attrFilter.setDefaultValue(precice::logging::BackendConfiguration::default_filter);
  tagSink.addAttribute(attrFilter);

  tagLog.addSubtag(tagSink);
  parent.addSubtag(tagLog);
}

void LogConfiguration::xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace("xmlTagCallback()", tag.getFullName());
  if (tag.getName() == "sink") {
    precice::logging::BackendConfiguration config;
    config.setOption("type", tag.getStringAttributeValue("type"));
    config.setOption("output", tag.getStringAttributeValue("output"));
    config.setOption("filter", tag.getStringAttributeValue("filter"));
    config.setOption("format", tag.getStringAttributeValue("format"));
    _logconfig.push_back(config);
  }
}

void LogConfiguration::xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace("xmlEndTagCallback()", tag.getFullName());
  if (tag.getName() == "log")
    precice::logging::setupLogging(_logconfig);
}

}} // namespace precice, config
