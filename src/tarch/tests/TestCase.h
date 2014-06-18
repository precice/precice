// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_TESTS_TESTCASE_H_
#define _TARCH_TESTS_TESTCASE_H_

#include "tarch/tests/TestMacros.h"
#include "tarch/compiler/CompilerSpecificSettings.h"
#include <string>


namespace tarch {
  namespace tests {
    class TestCase;
  }
}


/**
 * Represents one test case. Every test case should be a subclass of this
 * class implementing run(). Furthermore subtypes should use the assertion
 * macros included to check any assumption. All the test cases are managed by the
 * TestCaseCollection.
 *
 * Note that there is the subclass TestCaseWithScenario that has additional
 * creation methods for several (fluid) scenarios as well as some
 * comparison methods.
 *
 * If you want setup test case data before you run your tests, you have to
 * overwrite setUp(). To clean up your test fixture use tearDown(). Both
 * operations are called before run() is invoked. run() on the other hand is
 * called if and only if setUp() has been successful. To overwrite both fixture
 * management routines is optional - if the operations are not redefined they
 * equal nop.
 *
 * \see TestCaseCollection
 * \see TestCaseWithScenario
 *
 * !!! Rationale
 *
 * Whenever you write a test case implementation, we recommend to add the
 * code block
 * \code
#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",off)
#endif
   \endcode
 * at the very beginning of your implementation and
 * \code
#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
   \endcode
 * at the end of the implementation file.
 * @author Tobias Weinzierl, Wolfgang Eckhardt
 * @version $Revision: 1.31 $
 *
 * @todo Doku fuer macros vervollstaendigen!
 */
class tarch::tests::TestCase {
  protected:
    /*
     * this flag is used for the run_macro to indicate if an error has occurred
     * during the execution of a test method
     */
    bool _error;

    /**
     * Name of the test case.
     */
    std::string _testCaseName;

    /**
     * Error counter. This counter is increased on errors iff you use the
     * validate macros within your test cases.
     */
    int _errors;

    static std::string _outputDirectory;

    TestCase();

  public:
    /**
     * Constructor.
     *
     * @param testCaseName Name of the test case. If a class is tested, this
     *                     string should equal the class name added the
     *                     namespace.
     */
    TestCase( const std::string& testCaseName );

    /**
     * Destructor.
     */
    virtual ~TestCase();

    /**
     * @return Number of errors.
     */
    int getNumberOfErrors() const;

    /**
     * @return the name of the testcase
     */
    const std::string& getTestCaseName() const;

    void setTestCaseName(const std::string& testCaseName );

    /**
     * This routine is triggered by the TestCaseCollection
     */
    virtual void run() = 0;

    /**
     * Setup your test case.
     */
    void virtual setUp() = 0;

    /**
     * Set the output directory.
     */
    static void setOutputDirectory(const std::string & outputDirectory);
};

#endif
