// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_TESTS_TESTCASEREGISTRY_H_
#define _TARCH_TESTS_TESTCASEREGISTRY_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include <string>
#include "tarch/tests/TreeTestCaseCollection.h"

namespace tarch {
  namespace tests {
    class TestCaseRegistry;
    class TestCase;
  }
}


/**
 * TestCaseFactories are registered themselves with this class. It holds all
 * the tests, i.e. manages the test case collection tree.
 *
 * @author Wolfgang Eckhardt
 * @author Tobias Weinzierl
 */
class tarch::tests::TestCaseRegistry {
  private:
    TestCaseRegistry();

    TreeTestCaseCollection _globalTestCase;
    TreeTestCaseCollection _globalIntegrationTestCase;
  public:
    virtual ~TestCaseRegistry();

    static TestCaseRegistry& getInstance();

    void addTestCase(const std::string& testCaseName, TestCase* testCases);
    void addIntegrationTestCase(const std::string& testCaseName, TestCase* testCases);

    TestCase& getTestCaseCollection();
    TestCase& getIntegrationTestCaseCollection();
};

#endif /* TESTCASEREGISTRY_H_ */
