#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/tests/configurations/IntegrationTestConfiguration.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "utils/assertion.hpp"


#include "tarch/configuration/TopLevelConfigurationFactory.h"
registerTopLevelConfiguration(tarch::tests::configurations::IntegrationTestConfiguration)



precice::logging::Logger tarch::tests::configurations::IntegrationTestConfiguration::_log("tarch::tests::configurations::IntegrationTestConfiguration");


tarch::tests::configurations::IntegrationTestConfiguration::IntegrationTestConfiguration():
  _isValid(false),
  _outputDirectoryForTempFiles("") {
}


tarch::tests::configurations::IntegrationTestConfiguration::~IntegrationTestConfiguration() {
}


std::string tarch::tests::configurations::IntegrationTestConfiguration::getTag() const {
  return "run-integration-tests";
}


void tarch::tests::configurations::IntegrationTestConfiguration::toXML(std::ostream& out) const {
  out << "<!--" << std::endl
      << "  This is the configuration tag corresponding to tarch::tests::configurations::IntegrationTestConfiguration. " << std::endl
      << "  The xml tag has to contain the following attributes: " << std::endl
      << std::endl
      << "    | attribute name   | semantics | allowed values " << std::endl
      << "    | output-directory | Integration tests might write files to some place. Please specify where they are allowed to do so. | Any valid directory." << std::endl
      << std::endl
      << "  The xml tag may hold the following subtags:"
      << "    - " << _logConfiguration.getTag() << std::endl
      << "    - " << _logFormatConfiguration.getTag() << std::endl
      << "  -->" << std::endl;
  if ( isValid() ) {
    out << "<" << getTag() << " output-directory=\"" << _outputDirectoryForTempFiles << "\" >" << std::endl;
  }
  else {
    out << "<" << getTag() << " output-directory=\"a path\" >" << std::endl;
  }
 _logConfiguration.toXML(out);
 _logFormatConfiguration.toXML(out);
  out << "</" << getTag() << ">" << std::endl;
}


bool tarch::tests::configurations::IntegrationTestConfiguration::isValid() const {
	
	bool ret = _isValid
      && _logConfiguration.isValid()
      && _logFormatConfiguration.isValid();
	
  return ret;
}

void tarch::tests::configurations::IntegrationTestConfiguration::parseSubtag( precice::xml::Parser::CTag *pTag)
{	
	_isValid = true;

	if ( pTag->m_aAttributes.find("output-directory") == pTag->m_aAttributes.end() ) {
	_isValid = false;
	WARN("Missing or invalid attribute \"output-directory\" for tag <" + getTag() + ">");
	}
	else {
	_outputDirectoryForTempFiles = pTag->m_aAttributes["output-directory"];
	}


	for(auto tag : pTag->m_aSubTags)
	{
		if ( tag->m_Name == _logConfiguration.getTag() ) {
			_logConfiguration.parseSubtag(tag);
		}
		if ( tag->m_Name == _logFormatConfiguration.getTag() ) {
			_logFormatConfiguration.parseSubtag(tag);
		}
	}

	if (!_logConfiguration.isValid()) {
	WARN("subtag <" + _logConfiguration.getTag() + "> missing or invalid." );
	}
	if (!_logFormatConfiguration.isValid()) {
	WARN("subtag <" + _logFormatConfiguration.getTag() + "> missing or invalid." );
	}
}


tarch::configuration::TopLevelConfiguration* tarch::tests::configurations::IntegrationTestConfiguration::clone() const {
  return new IntegrationTestConfiguration();
}


int tarch::tests::configurations::IntegrationTestConfiguration::interpreteConfiguration() {
  TestCase::setOutputDirectory( _outputDirectoryForTempFiles );
  TestCaseRegistry::getInstance().getIntegrationTestCaseCollection().run();
  return TestCaseRegistry::getInstance().getIntegrationTestCaseCollection().getNumberOfErrors();
}
