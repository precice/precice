// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Configuration.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace config {

tarch::logging::Log Configuration:: _log("precice::config::Configuration");

//const std::string& Configuration:: getTag()
//{
//  static std::string tag("precice-configuration");
//  return tag;
//}

Configuration:: Configuration()
:
  _tag(*this, "precice-configuration", utils::XMLTag::OCCUR_ONCE),
  //_isValid(false),
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

//  using namespace utils;
  //XMLTag tag(getTag(), XMLTag::OCCUR_ONCE);
  //XMLTag tagLogFilter(*this, _logFilterConfig.getTag(), XMLTag::OCCUR_ARBITRARY);
  //XMLTag tagLogFormat(*this, _logFormatConfig.getTag(), XMLTag::OCCUR_NOT_OR_ONCE);
  //XMLTag tagSolverInterface(*this_solverInterfaceConfig.getTag(), XMLTag::OCCUR_ONCE);
  //_tag.addSubtag(tagLogFilter);
  //_tag.addSubtag(tagLogFormat);
  //_tag.addSubtag(tagSolverInterface);
}

utils::XMLTag& Configuration:: getXMLTag()
{
  return _tag;
}

//bool Configuration:: parseSubtag
//(
//  utils::XMLTag::XMLReader* xmlReader )
//{
//  _isValid = _tag.parse(xmlReader);
//  return _isValid;
//}

void Configuration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlTagCallback()", tag.getName());
//  if (tag.getName() == getTag()){
//    // Nothing to do
//  }
//  else if (callingTag.getName() == _logFilterConfig.getTag()){
//    _logFilterConfig.parseSubtag(xmlReader);
//    return _logFilterConfig.isValid();
//  }
//  else if (callingTag.getName() == _logFormatConfig.getTag()){
//    _logFormatConfig.parseSubtag(xmlReader);
//    return _logFormatConfig.isValid();
//  }
//  else {
//    preciceError("xmlTagCallback()", "Received callback from tag " << callingTag.getName());
//  }
//  else if (callingTag.getName() == _solverInterfaceConfig.getTag()){
//    return _solverInterfaceConfig.parseSubtag(xmlReader);
//  }
}

void Configuration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
}

//bool Configuration:: isValid() const
//{
//  return _isValid;
//}

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
