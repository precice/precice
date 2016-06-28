#include "LogConfiguration.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"

namespace precice {
namespace config {

/// Logging device.
precice::logging::Logger precice::config::LogConfiguration::_log("logging::config::LogConfiguration");

LogConfiguration::LogConfiguration
(
  utils::XMLTag& parent)
:
  TAG("log"),
  ATTR_FILE("file")
{
  using namespace utils;
  std::string doc;
  XMLTag tag(*this, TAG, XMLTag::OCCUR_ONCE);
  doc = "Configures logging";
  tag.setDocumentation(doc);

  XMLAttribute<std::string> attrFile(ATTR_FILE);
  doc = "Logging configuration file.";
  attrFile.setDocumentation(doc);
  attrFile.setDefaultValue("log.conf");

  tag.addAttribute(attrFile);

  parent.addSubtag(tag);
}

void LogConfiguration::xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace("xmlTagCallback()", tag.getFullName());
  if (tag.getName() == TAG)
    precice::logging::setupLogging(tag.getStringAttributeValue(ATTR_FILE));
  else
    precice::logging::setupLogging(); // log tag has not been supplied., use default log config file
}

void LogConfiguration::xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace("xmlEndTagCallback()", tag.getFullName());
}

}} // namespace precice, config
