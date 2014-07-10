#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#include "tarch/la/tests/ScalarTest.h"

#include "tarch/la/ScalarOperations.h"
#include "tarch/tests/TestCaseFactory.h"

#include "tarch/compiler/CompilerSpecificSettings.h"


registerTest(tarch::la::ScalarTest)

namespace tarch {
namespace la {

ScalarTest::ScalarTest ()
:
  TestCase("tarch::la::ScalarTest")
{}

void ScalarTest::run ()
{
  testMethod (testComparison);
  testMethod (testAbs);
}

void ScalarTest::testComparison()
{
  // @todo Whoever wrote this has never tested it with CLX
  #if !defined(CompilerCLX)
  double a = 1.0;
  double b = 2.0;
  double eps = 1e-14;

  validate (greater(b, a, eps));
  validate (not greater(a, a - eps, eps));
  validate (greater(a, a - 10.0 * eps, eps));

  validate (not greaterEquals(a, b, eps));
  validate (greaterEquals(b, a, eps));
  validate (greaterEquals(a, a, eps));
  validate (greaterEquals(a, a + 0.1*eps, eps));
  validate (not greaterEquals(a, a + 10*eps, eps));

  validate (smaller(a, b, eps));
  validate (not smaller(a, a + eps, eps));
  validate (smaller(a, a + 10.0 * eps, eps));

  validate (smallerEquals(a, b, eps));
  validate (smallerEquals(a, a + eps, eps));
  validate (smallerEquals(a + eps, a, eps));
  validate (smallerEquals(a, a + 10.0 * eps, eps));

  validate (not equals(a, b, eps));
  validate (equals(a, a, eps));
  validate (equals(a, a + eps, eps));
  validate (not equals(a, a + 10.0 * eps, eps));
  validate (equals(a, a + 10.0 * eps, 10.0 * eps));
  #endif
}

void ScalarTest::testAbs()
{
  double a = -1.0;
  int b = -1;

  validateEquals (abs(a), 1.0);
  validateEquals (abs(b), 1);
}

}} // namespace tarch, la
