// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_TESTS_TREETESTCASECOLLECTION_H_
#define _TARCH_TESTS_TREETESTCASECOLLECTION_H_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

#include <list>
#include <map>

namespace tarch {
  namespace tests {
    class TreeTestCaseCollection;
  }
}


/**
 *
 * @author Tobias Weinzierl
 * @author Wolfgang Eckhardt
 * @version $Revision: 1.7 $
 */
class tarch::tests::TreeTestCaseCollection: public tarch::tests::TestCase {
  private:
    /**
     * Sequence of test cases that are executes on a run() call. The class is
     * not responsible for destroying them.
     */
    std::list<TestCase*>                            _testCases;
    std::map<std::string, TreeTestCaseCollection*>  _subTests;

    /**
     * Tells the collection wheather to log a progression bar.
     */
    bool _writeToLog;

    /**
     * Log interface the class writes to.
     */
    static precice::logging::Logger _log;

    /**
     * Tells whether the testcases contained should be deleted uppon destruction of this object.
     */
    bool _deleteTestCases;

    static bool isNameWithoutHierarchy(const std::string& testCaseName);
    static std::string getFirstIdentifierInHierarchy(const std::string& testCaseName);
    static std::string getRemainingPathWithoutIdentifier(const std::string& testCaseName);
  public:

    /**
     * Creates a test case collection.
     *
     * @param testCaseCollectionName Name of the test case collection. May not contain ::
     */
    TreeTestCaseCollection(
      const std::string& testCaseCollectionName,
      bool deleteTestCases,
      bool writeToLog
    );

    /**
     * Destructor
     */
    virtual ~TreeTestCaseCollection();

    /**
     * Runs all test cases assigned.
     */
    virtual void run();

    /**
     * Runs all test cases assigned.
     */
    virtual void run( const std::string& prefix );

    /**
     * Call setUp() on all associated sub test cases.
     */
    virtual void setUp();

    /**
     * Adds a new test case.
     *
     * Although you pass a pointer, you are still responsible for the instance
     * management.
     *
     * You have to specify in the constructor of a TestCaseCollection whether
     * it should delete testcases contained on construction.
     *
     * @Note: default is set to false (i.e. no destruction)
     */
    void addTestCase( const std::string& path, const std::string& fullQualifiedPath, TestCase* testCase );
};


#endif
