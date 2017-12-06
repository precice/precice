#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/configuration/ConfigurationRegistry.h"
#include "tarch/configuration/TopLevelConfiguration.h"
//#include "utils/assertion.hpp"
//#include "logging/LogMacros.hpp"
#include "utils/Globals.hpp"

precice::logging::Logger tarch::configuration::ConfigurationRegistry::_log("tarch::configuration::ConfigurationRegistry");

tarch::configuration::ConfigurationRegistry::ConfigurationRegistry():
//_log("tarch::configuration::ConfigurationRegistry"),
_topLevelTags(){
}

std::list<tarch::configuration::TopLevelConfiguration*> 
tarch::configuration::ConfigurationRegistry::parseTag(const std::string& filename, const std::string& topLevelTag)
{
	std::list<tarch::configuration::TopLevelConfiguration*>  result;
  
	precice::xml::Parser p(filename);
	parseTag(p.getRootTag(), topLevelTag, result);
	
	return result;
}

void tarch::configuration::ConfigurationRegistry::parseTag(precice::xml::Parser::CTag *pTag, const std::string& topLevelTag,
std::list<tarch::configuration::TopLevelConfiguration*> &result)
{
	bool error = false;
	
	if(!pTag->m_aSubTags.empty())
	{
		std::string currentTag = pTag->m_Name;
		
		if (topLevelTag == currentTag) {
		}
		else if (_topLevelTags.find(currentTag)==_topLevelTags.end()) 
		{
			WARN("invalid top level tag. Received <" + currentTag + "> but expected <"
				+ topLevelTag + "> or any top level tag"
			);

			enlistAvailableTopLevelTags();
			error = true;
		}
		else 
		{
			TopLevelConfiguration* configuration = _topLevelTags[currentTag]->clone();

			configuration->parseSubtag(pTag);

			if (configuration->isValid()) {
			  result.push_back(configuration);
			}
			else {
				WARN("Invalid subtag: " + pTag->m_Name);
				error = true;
			}
		}
	}
	
	if(error)
		ERROR("invalid tags were used");
	
	
	if(!pTag->m_aSubTags.empty())
	{
		for(auto subtags : pTag->m_aSubTags)
		{
			precice::xml::Parser::CTag *pSubTag = (precice::xml::Parser::CTag *)subtags;
			
			parseTag(pSubTag, topLevelTag, result);
		}
	}
}


tarch::configuration::ConfigurationRegistry::~ConfigurationRegistry() {
}


tarch::configuration::ConfigurationRegistry& tarch::configuration::ConfigurationRegistry::getInstance() {
  static ConfigurationRegistry singleton;
  return singleton;
}


void tarch::configuration::ConfigurationRegistry::addTopLevelConfiguration( TopLevelConfiguration* configuration ) {
  assertion( _topLevelTags.count(configuration->getTag())==0, configuration->getTag() );
  _topLevelTags[ configuration->getTag() ] = configuration;
}

void tarch::configuration::ConfigurationRegistry::addTopLevelConfigurationFactory( BaseTopLevelConfigurationFactory* factory ) {
  _topLevelTagFactories.push_back(factory);
}

void tarch::configuration::ConfigurationRegistry::initTopLevelConfigurationFactories() {
  for(auto & elem : _topLevelTagFactories) {
    elem->init();
  }
}

void tarch::configuration::ConfigurationRegistry::writeDummyConfigFile(std::ostream& out) const {
  std::cout << "<your-top-level-tag>" << std::endl
            << "<!--" << std::endl
            << "  This file enlists all the tags that are available. " << std::endl
            << "  ->" << std::endl;
  for (const auto & elem : _topLevelTags) {
    elem.second->toXML(out);
  }
  std::cout << "</your-top-level-tag>" << std::endl;
}


void tarch::configuration::ConfigurationRegistry::enlistAvailableTopLevelTags() const {
  std::ostringstream msg;
  msg << "available top level tags: ";
  for (const auto & elem : _topLevelTags) {
    msg << "<" << elem.first << "> ";
  }
  INFO(msg.str() );
}


void tarch::configuration::ConfigurationRegistry::freeConfigurations(std::list<TopLevelConfiguration*>& configurations) {
  for (auto & configuration : configurations) {
    delete configuration;
  }
}


std::list<tarch::configuration::TopLevelConfiguration*>
tarch::configuration::ConfigurationRegistry::readFile(
  const std::string& filename,
  const std::string& topLevelTag
) {  
  return parseTag(filename, topLevelTag);
}
