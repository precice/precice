// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_TESTS_TEST_MACROS_H_
#define _TARCH_TESTS_TEST_MACROS_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include <cmath>
#include <iostream>

/**
 * Run a test method and check for errors. It is necessary to
 * use assertion-macros which set member _error in a TestCase.
 *
 * @param name the name of the method to call (without parentheses!)
 */
#define testMethod(name) { \
      _error = false; \
      name(); \
      if (_error) { \
        std::cerr << "testcase " << _testCaseName << "::" << #name << "() failed" << std::endl << std::endl; \
      } \
    };

/**
 * Automatically register a testCase.
 */
#define registerTest(name) \
  static tarch::tests::TestCaseFactory<name> thisTestCaseFactoryInstance( tarch::tests::TestCaseFactory<name>::UnitTest, #name );


/**
 * Automatically register an IntegrationTest.
 */
#define registerIntegrationTest(name) \
  static tarch::tests::TestCaseFactory<name> thisTestCaseFactoryInstance( tarch::tests::TestCaseFactory<name>::IntegrationTest, #name );


#define testNumericalEquals(lhs,rhs) \
  (std::abs((rhs) -(lhs)) <= 1.0e-10)


#define validate(booleanExpr) if (!(booleanExpr)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  boolean test failed " << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #booleanExpr << std::endl; \
    return;\
  }

/**
 * parameter message has to be a std::string or char*
 */
#define validateWithMessage(booleanExpr,message) if (!(booleanExpr)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  boolean test failed " << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #booleanExpr << std::endl \
              << "  message: " << message << std::endl; \
    return;\
  }

#define validateWithParams1(booleanExpr,param0) if (!(booleanExpr)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  boolean test failed " << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #booleanExpr << std::endl \
              << "  parameter " << #param0 << "=" << param0 << std::endl; \
    return;\
  }

#define validateWithParams2(booleanExpr,param0,param1) if (!(booleanExpr)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  boolean test failed " << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #booleanExpr << std::endl \
              << "  parameter " << #param0 << "=" << param0 << std::endl \
              << "  parameter " << #param1 << "=" << param1 << std::endl; \
    return;\
  }

#define validateWithParams3(booleanExpr,param0,param1,param2) if (!(booleanExpr)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  boolean test failed " << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #booleanExpr << std::endl \
              << "  parameter " << #param0 << "=" << param0 << std::endl \
              << "  parameter " << #param1 << "=" << param1 << std::endl \
              << "  parameter " << #param2 << "=" << param2 << std::endl; \
    return;\
  }

#define validateWithParams4(booleanExpr,param0,param1,param2,param3) if (!(booleanExpr)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  boolean test failed " << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #booleanExpr << std::endl \
              << "  parameter " << #param0 << "=" << param0 << std::endl \
              << "  parameter " << #param1 << "=" << param1 << std::endl \
              << "  parameter " << #param2 << "=" << param2 << std::endl \
              << "  parameter " << #param3 << "=" << param3 << std::endl; \
    return;\
  }

#define validateEquals(actualValue, validValue) if (!(actualValue == validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl; \
    return;\
  }

#define validateEqualsWithMessage(actualValue, validValue,message) if (!(actualValue == validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  message " <<  message << std::endl; \
    return;\
  }

#define validateEqualsWithParams1(actualValue, validValue, param0) if (!(actualValue == validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << std::endl; \
    return;\
  }

#define validateEqualsWithParams2(actualValue, validValue, param0, param1) if (!(actualValue == validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << std::endl; \
    return;\
  }


#define validateEqualsWithParams3(actualValue, validValue, param0, param1, param2) if (!(actualValue == validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr  << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << std::endl; \
    return;\
  }

#define validateEqualsWithParams4(actualValue, validValue, param0, param1, param2, param3) if (!(actualValue == validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << ", parameter " << #param3 << "=" << param3 << std::endl; \
    return;\
  }


#define validateNotEqual(actualValue, validValue) if ((actualValue == validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  inequality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "!="  << #validValue << std::endl; \
    return;\
  }

#define validateNumericalEquals(actualValue, validValue) if (!testNumericalEquals(actualValue,validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
    return;\
  }

#define validateNumericalEqualsWithEps(actualValue, validValue, eps) \
  if (  (fabs((actualValue) - (validValue)) > eps ) ) { \
    _errors++; \
    _error = true; \
    std::cerr << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
    return;\
  }

#define validateNumericalEqualsWithEpsWithParams1(actualValue, validValue, eps, param0) \
  if (  (fabs((actualValue) - (validValue)) > eps ) ) { \
    _errors++; \
    _error = true; \
    std::cerr<< "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
    return;\
  }

#define validateNumericalEqualsWithParams1(actualValue, validValue, param0) if (!testNumericalEquals(actualValue,validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
    return;\
  }

#define validateNumericalEqualsWithParams2(actualValue, validValue, param0, param1) if (!testNumericalEquals(actualValue,validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
    return;\
  }

#define validateNumericalEqualsWithParams3(actualValue, validValue, param0, param1, param2) if (!testNumericalEquals(actualValue,validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
    return;\
  }

#define validateNumericalEqualsWithParams4(actualValue, validValue, param0, param1, param2, param3) if (!testNumericalEquals(actualValue,validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << ", parameter " << #param3 << "=" << param3 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
    return;\
  }

#define validateNumericalEqualsWithParams5(actualValue, validValue, param0, param1, param2, param3, param4) if (!testNumericalEquals(actualValue,validValue)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << ", parameter " << #param3 << "=" << param3<< ", parameter " << #param4 << "=" << param4 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
    return;\
  }

/**
 * Macro for validating the (numerical) equality of two Vectors.
 *
 * The "difference" is here implemented as Euclidean norm of the difference
 * vector.
 *
 * @param actual             Vector under test.
 * @param valid              Valid reference vector.
 * @param testCaseMethodName String containing the name of the test case method.
 */
#define validateNumericalVectorEquals(actual, valid) \
  if (!tarch::la::equals(actual, valid)) { \
    _errors++; \
    _error = true; \
    std::cerr << "  numerical vector equality test failed: " << actual << " instead of " << valid << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actual << "=="  << #valid << std::endl; \
    return;\
  }


#endif /* VALIDATEMACROS_H_ */
