#include "Configuration.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace config {

tarch::logging::Log Configuration:: _log("precice::config::Configuration");

Configuration:: Configuration()
:
  _tag(*this, "precice-configuration", utils::XMLTag::OCCUR_ONCE),
  _logFilterConfig(_tag),
  _logFormatConfig(_tag),
  _solverInterfaceConfig(_tag)
{
  _tag.setDocumentation("Main tag containing preCICE configuration.");
  _tag.addNamespace("data");
  _tag.addNamespace("geometry");
  _tag.addNamespace("spacetree");
  _tag.addNamespace("communication");
  _tag.addNamespace("mapping");
  _tag.addNamespace("export");
  _tag.addNamespace("action");
  _tag.addNamespace("server");
  _tag.addNamespace("coupling-scheme");
  _tag.addNamespace("post-processing");
}

utils::XMLTag& Configuration:: getXMLTag()
{
  return _tag;
}

void Configuration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlTagCallback()", tag.getName());
}

void Configuration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlEndTagCallback()", tag.getName());
}

const LogFilterConfiguration& Configuration:: getLogFilterConfiguration() const
{
  return _logFilterConfig;
}

const LogOutputFormatConfiguration& Configuration:: getLogFormatConfiguration() const
{
  return _logFormatConfig;
}

const SolverInterfaceConfiguration&
Configuration:: getSolverInterfaceConfiguration() const
{
  return _solverInterfaceConfig;
}

}} // namespace precice, config
