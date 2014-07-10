// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_TESTS_TESTCASECOLLECTION_H_
#define _TARCH_TESTS_TESTCASECOLLECTION_H_

#ifdef Parallel
#include <mpi.h>
#endif

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

#include <list>

namespace tarch {
  namespace tests {
    class TestCaseCollection;
  }
}


/**
 * Contains a sequence of tests which have to be executed sequentially. The
 * object manages a sequence of references. The user is responsible for
 * creating and destroying the test case instances.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.7 $
 */
class tarch::tests::TestCaseCollection: public tarch::tests::TestCase {
  private:
    /**
     * Sequence of test cases that are executes on a run() call. The class is
     * not responsible for destroying them.
     */
    std::list<TestCase*>            _testCases;

    /**
     * Tells the collection wheather to log a progression bar.
     */
    bool _writeToLog;

    /**
     * Log interface the class writes to.
     */
    static tarch::logging::Log _log;

    /**
     * Tells whether the testcases contained should be deleted uppon destruction of this object.
     */
    bool _deleteTestCases;

    TestCaseCollection();
  public:
    /**
     * Creates a test case collection.
     *
     * @param testCaseCollectionName Name of the test case collection.
     */
    TestCaseCollection(const std::string& testCaseCollectionName, bool deleteTestCases = false, bool writeToLog = true);

    /**
     * Destructor
     */
    virtual ~TestCaseCollection();

    /**
     * Runs all test cases assigned.
     */
    virtual void run();

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
    void addTestCase( TestCase* testCase );
};


#endif
