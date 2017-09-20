#include "Configuration.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace config {

logging::Logger Configuration:: _log("config::Configuration");

Configuration:: Configuration()
:
  _tag(*this, "precice-configuration", utils::XMLTag::OCCUR_ONCE),
  _logConfig(_tag),
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
  TRACE(tag.getName());
}

void Configuration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  TRACE(tag.getName());
}

const SolverInterfaceConfiguration&
Configuration:: getSolverInterfaceConfiguration() const
{
  return _solverInterfaceConfig;
}

}} // namespace precice, config
