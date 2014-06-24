#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/tests/TestCase.h"


std::string tarch::tests::TestCase::_outputDirectory = "./";

tarch::tests::TestCase::TestCase()
:
  _testCaseName(),
  _errors(0)
{}

void tarch::tests::TestCase::setTestCaseName(const std::string& name) {
  _testCaseName = name;
}

void tarch::tests::TestCase::setOutputDirectory(const std::string & outputDirectory) {
  _outputDirectory = outputDirectory;
  if ( _outputDirectory.at( _outputDirectory.length()-1 )!='/' ) {
    _outputDirectory.append( "/" );
  }
}


tarch::tests::TestCase::TestCase( const std::string& testCaseName ):
  _testCaseName(testCaseName),
  _errors(0) {
}

tarch::tests::TestCase::~TestCase() {}

int tarch::tests::TestCase::getNumberOfErrors() const {
  return _errors;
}

const std::string& tarch::tests::TestCase::getTestCaseName() const {
  return _testCaseName;
}
