#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/tests/configurations/IntegrationTestConfiguration.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/Assertions.h"


#include "tarch/configuration/TopLevelConfigurationFactory.h"
registerTopLevelConfiguration(tarch::tests::configurations::IntegrationTestConfiguration)



tarch::logging::Log tarch::tests::configurations::IntegrationTestConfiguration::_log("tarch::tests::configurations::IntegrationTestConfiguration");


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
  return _isValid
      && _logConfiguration.isValid()
      && _logFormatConfiguration.isValid();
}


void tarch::tests::configurations::IntegrationTestConfiguration::parseSubtag( tarch::irr::io::IrrXMLReader* xmlReader ) {
  assertion( xmlReader != 0 );

  _isValid = true;

  if ( xmlReader->getAttributeValue("output-directory")==0 ) {
    _isValid = false;
    _log.error(
      "parseSubtag(...)",
      "missing or invalid attribute \"output-directory\" for tag <" + getTag() +
      ">"
    );
  }
  else {
    _outputDirectoryForTempFiles = xmlReader->getAttributeValue("output-directory");
  }


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


tarch::configuration::TopLevelConfiguration* tarch::tests::configurations::IntegrationTestConfiguration::clone() const {
  return new IntegrationTestConfiguration();
}


int tarch::tests::configurations::IntegrationTestConfiguration::interpreteConfiguration() {
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntries( _logConfiguration.getFilterList() );
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( _logFormatConfiguration );
  TestCase::setOutputDirectory( _outputDirectoryForTempFiles );
  TestCaseRegistry::getInstance().getIntegrationTestCaseCollection().run();
  return TestCaseRegistry::getInstance().getIntegrationTestCaseCollection().getNumberOfErrors();
}
