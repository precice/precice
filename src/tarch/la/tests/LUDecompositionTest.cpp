#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#include "LUDecompositionTest.h"
#include "tarch/la/LUDecomposition.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/Matrix.h"
#include "tarch/la/DynamicColumnMatrix.h"
#include "tarch/la/Vector.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(tarch::la::LUDecompositionTest)

namespace tarch {
namespace la {

LUDecompositionTest::LUDecompositionTest()
:
  TestCase ("tarch::la::LUDecompositionTest")
{}

void LUDecompositionTest::run()
{
  testMethod (testLUNoPivoting);
  testMethod (testLU);
}

void LUDecompositionTest::testLUNoPivoting()
{
  Vector<3,int> pivots3(0);
  DynamicMatrix<double> dynA (3,3);
  assignList(dynA) =
      5.0, 2.0, 3.0,
      2.0, 3.0, 6.0,
      1.0, 2.0, 4.0;
  lu (dynA, pivots3);
  Matrix<3,3,double> result3by3;
  assignList(result3by3) =
      5.0, 2.0, 3.0,
      0.4, 2.2, 4.8,
      0.2, 0.7272727272727272, -0.0909090909090909;
  //TODO validateWithParams1 (equals(dynA, result3by3), dynA);
  Vector<3,int> validPivots3(0,1,2);
  validateEqualsWithParams1 (pivots3, validPivots3, pivots3);

  Vector<4,int> pivots4(0);
  DynamicColumnMatrix<double> dyncolA (4,4);
  assignList(dyncolA) =
      5.0, 2.0, 3.0, 1.0,
      2.0, 3.0, 6.0, 1.0,
      2.0, 2.0, 5.0, 4.0,
      1.0, 1.0, 1.0, 1.0;
  lu (dyncolA, pivots4);
  Matrix<4,4,double> result4by4;
  assignList(result4by4) =
      5.0, 2.0, 3.0, 1.0,
      0.4, 2.2, 4.8, 0.6,
      0.4, 0.5454545454545454, 1.1818181818181818, 3.2727272727272727,
      0.2, 0.2727272727272727, -0.7692307692307689, 3.153846153846153;
  Vector<4,int> validPivots4(0,1,2,3);
  //TODO validateWithParams1 (equals(dyncolA, result4by4), dyncolA);
  validateEqualsWithParams1 (pivots4, validPivots4, pivots4);
}

void LUDecompositionTest::testLU()
{
  Matrix<4,4,double> A;
  assignList(A) =
      5.0, 2.0, 3.0, 1.0,
      6.0, 3.0, 6.0, 1.0,
      3.0, 2.0, 5.0, 4.0,
      4.0, 1.0, 2.0, 1.0;
  Matrix<4,4,double> Acopy(A);
  Vector<4,int> pivots;
  lu(A,pivots);
  Matrix<4,4,double> L(A);
  for (int i=0; i < 4; i++){
    L(i,i) = 1.0;
    for (int j=i+1; j < 4; j++){
      L(i,j) = 0.0;
    }
  }
  Matrix<4,4,double> R(A);
  for (int j=0; j < 4; j++){
    for (int i=j+1; i < 4; i++){
      R(i,j) = 0.0;
    }
  }
  assign(A) = 0.0;
  multiply (L, R, A);
  for (int i=0; i < 4; i++){
    for (int j=0; j < 4; j++){
      double temp = Acopy(i,j);
      Acopy(i,j) = Acopy(pivots(i),j);
      Acopy(pivots(i),j) = temp;
    }
  }
  //TODO validateWithParams1 (equals(A, Acopy), A);
}

}} // namespace tarch, la
