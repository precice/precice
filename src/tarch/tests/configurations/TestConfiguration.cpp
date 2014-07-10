#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/tests/configurations/TestConfiguration.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/Assertions.h"


#include "tarch/configuration/TopLevelConfigurationFactory.h"
registerTopLevelConfiguration(tarch::tests::configurations::TestConfiguration)


tarch::logging::Log tarch::tests::configurations::TestConfiguration::_log("tarch::tests::configurations::TestConfiguration");


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
  return _isValid
      && _logConfiguration.isValid()
      && _logFormatConfiguration.isValid();
}


void tarch::tests::configurations::TestConfiguration::parseSubtag( tarch::irr::io::IrrXMLReader* xmlReader ) {
  assertion( xmlReader != 0 );

  _isValid = true;

  while(
    (xmlReader->getNodeType()!=irr::io::EXN_ELEMENT_END) &&
    (xmlReader->read() )
  ) {
    if ( xmlReader->getNodeType()==irr::io::EXN_ELEMENT ) {
      if ( xmlReader->getNodeName() == _logConfiguration.getTag() ) {
        _logConfiguration.parseSubtag(xmlReader);
      }
      if ( xmlReader->getNodeName() == _logFormatConfiguration.getTag() ) {
        _logFormatConfiguration.parseSubtag(xmlReader);
      }
    }
  }

  if (
    (xmlReader->getNodeType()!=irr::io::EXN_ELEMENT_END) ||
    (xmlReader->getNodeName()!=getTag())
  ) {
    _log.error(
      "parseSubtag(...)",
      "expected closing tag for " + getTag() +
      ", but received tag <" + xmlReader->getNodeName() + ">"
    );
    _isValid = false;
  }

  if (!_logConfiguration.isValid()) {
    _log.error( "parse(...)", "subtag <" + _logConfiguration.getTag() + "> missing or invalid." );
  }
  if (!_logFormatConfiguration.isValid()) {
    _log.error( "parse(...)", "subtag <" + _logFormatConfiguration.getTag() + "> missing or invalid." );
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
