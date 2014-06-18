#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/tests/TestCaseRegistry.h"

#include "tarch/tests/TestCase.h"


tarch::tests::TestCaseRegistry::TestCaseRegistry():
  _globalTestCase("", false, true),
  _globalIntegrationTestCase("", false, true) {
}


tarch::tests::TestCaseRegistry::~TestCaseRegistry() {
}


tarch::tests::TestCaseRegistry& tarch::tests::TestCaseRegistry::getInstance() {
  static tarch::tests::TestCaseRegistry singleton;
  return singleton;
}


void tarch::tests::TestCaseRegistry::addTestCase(const std::string& testCaseName, tarch::tests::TestCase* testCases) {
  _globalTestCase.addTestCase( testCaseName, testCaseName, testCases );
}


void tarch::tests::TestCaseRegistry::addIntegrationTestCase(const std::string& testCaseName, tarch::tests::TestCase* testCases) {
  _globalIntegrationTestCase.addTestCase( testCaseName, testCaseName, testCases );
}


tarch::tests::TestCase& tarch::tests::TestCaseRegistry::getTestCaseCollection() {
  _globalTestCase.setUp();
  return _globalTestCase;
}


tarch::tests::TestCase& tarch::tests::TestCaseRegistry::getIntegrationTestCaseCollection() {
  _globalIntegrationTestCase.setUp();
  return _globalIntegrationTestCase;
}
