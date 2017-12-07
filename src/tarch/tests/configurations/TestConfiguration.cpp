#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/tests/configurations/TestConfiguration.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "utils/assertion.hpp"


#include "tarch/configuration/TopLevelConfigurationFactory.h"
registerTopLevelConfiguration(tarch::tests::configurations::TestConfiguration)


precice::logging::Logger tarch::tests::configurations::TestConfiguration::_log("tarch::tests::configurations::TestConfiguration");


tarch::tests::configurations::TestConfiguration::TestConfiguration():
  tarch::configuration::TopLevelConfiguration(),
  _isValid(false) {
}


tarch::tests::configurations::TestConfiguration::~TestConfiguration() {
}


std::string tarch::tests::configurations::TestConfiguration::getTag() const {
  return "run-tests";
}


void tarch::tests::configurations::TestConfiguration::toXML(std::ostream& out) const {
  out << "<!--" << std::endl
      << "  This is the configuration tag corresponding to tarch::tests::configurations::IntegrationTestConfiguration. " << std::endl
      << "  The xml tag may not contain any attributes. " << std::endl
      << std::endl
      << std::endl
      << "  The xml tag may hold the following subtags:"
      << "    - " << _logConfiguration.getTag() << std::endl
      << "    - " << _logFormatConfiguration.getTag() << std::endl
      << "  -->" << std::endl;
  out << "<" << getTag() << " >" << std::endl;
 _logConfiguration.toXML(out);
 _logFormatConfiguration.toXML(out);
  out << "</" << getTag() << ">" << std::endl;
}


bool tarch::tests::configurations::TestConfiguration::isValid() const {
	bool ret = _isValid
      && _logConfiguration.isValid()
      && _logFormatConfiguration.isValid();
	   
  return ret;
}

void tarch::tests::configurations::TestConfiguration::parseSubtag(precice::xml::ConfigParser::CTag *pTag)
{
	_isValid = true;
	
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


tarch::configuration::TopLevelConfiguration* tarch::tests::configurations::TestConfiguration::clone() const {
  return new TestConfiguration();
}


int tarch::tests::configurations::TestConfiguration::interpreteConfiguration() {
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntries( _logConfiguration.getFilterList() );
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( _logFormatConfiguration );
  TestCaseRegistry::getInstance().getTestCaseCollection().run();
  return TestCaseRegistry::getInstance().getTestCaseCollection().getNumberOfErrors();
}
