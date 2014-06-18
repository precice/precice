#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#include "TraitsTest.h"
#include "tarch/la/traits/EqualScalars.h"
#include "tarch/la/Vector.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(tarch::la::TraitsTest)

namespace tarch {
namespace la {

TraitsTest::TraitsTest ()
:
  TestCase("tarch::la::TraitsTest")
{}

void TraitsTest::run()
{
  testMethod (testEqualScalars);
}

//void TraitsTest::testIsEqual()
//{
//  bool equal = IsEqual<int,int>::value;
//  bool notEqual = not IsEqual<int,double>::value;
//  validate (equal);
//  validate (notEqual);
//}

void TraitsTest:: testEqualScalars()
{
  bool equal = EqualScalars<Vector<3,double>,Vector<2,double> >::value;
  validate (equal);
  equal = EqualScalars<Vector<3,double>,DynamicVector<double> >::value;
  validate (equal);
  equal = EqualScalars<Vector<3,double>,Vector<3,int> >::value;
  validate (not equal);
  equal = EqualScalars<Vector<3,float>,DynamicVector<double> >::value;
  validate (not equal);
}

}} // namespace tarch, la
