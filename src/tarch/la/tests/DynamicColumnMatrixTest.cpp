#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#include "DynamicColumnMatrixTest.h"
#include "tarch/la/DynamicColumnMatrix.h"
#include <string>
#include <sstream>

#include "tarch/tests/TestCaseFactory.h"
registerTest(tarch::la::DynamicColumnMatrixTest)

namespace tarch {
namespace la {

DynamicColumnMatrixTest::DynamicColumnMatrixTest ()
:
  TestCase ("tarch::la::MatrixTest")
{}

void DynamicColumnMatrixTest::run ()
{
  testMethod (testBasics);
}

void DynamicColumnMatrixTest::testBasics ()
{
  DynamicColumnMatrix<int> matrix(2, 3, 5);
  validateEquals (matrix.size(), 6);
  validateEquals (matrix.rows(), 2);
  validateEquals (matrix.cols(), 3);
  validateEquals (matrix(0,0), 5);
  validateEquals (matrix(0,1), 5);
  validateEquals (matrix(1,2), 5);
  assignList(matrix) = 1, 2, 3,
                       4, 5, 6;
  validateEquals (matrix(0,0), 1);
  validateEquals (matrix(0,1), 2);
  validateEquals (matrix(0,2), 3);
  validateEquals (matrix(1,0), 4);
  validateEquals (matrix(1,1), 5);
  validateEquals (matrix(1,2), 6);
  validateEquals (matrix.column(0).size(), 2);
  validateEquals (matrix.column(1).size(), 2);
  validateEquals (matrix.column(2).size(), 2);
  validateEquals (matrix.column(0)[0], 1);
  validateEquals (matrix.column(1)[0], 2);
  validateEquals (matrix.column(2)[1], 6);

  DynamicColumnMatrix<int> matrix2(matrix);
  validateEquals (matrix2.size(), 6);
  validateEquals (matrix2.rows(), 2);
  validateEquals (matrix2.cols(), 3);
  validateEquals (matrix2(0,0), 1);
  validateEquals (matrix2(0,1), 2);
  validateEquals (matrix2(0,2), 3);
  validateEquals (matrix2(1,0), 4);
  validateEquals (matrix2(1,1), 5);
  validateEquals (matrix2(1,2), 6);
}

void DynamicColumnMatrixTest::testColumnManipulations ()
{
  DynamicColumnMatrix<int> matrix(2, 3);
  validateEquals (matrix.size(), 6);
  validateEquals (matrix.rows(), 2);
  validateEquals (matrix.cols(), 3);
  assignList(matrix) = 1, 2, 3,
                       4, 5, 6;

  // Remove columns
  matrix.remove(2);
  validateEquals (matrix.size(), 4);
  validateEquals (matrix.rows(), 2);
  validateEquals (matrix.cols(), 2);
  matrix.remove(1);
  validateEquals (matrix.size(), 2);
  validateEquals (matrix.rows(), 2);
  validateEquals (matrix.cols(), 1);
  matrix.remove(0);
  validateEquals (matrix.size(), 0);
  validateEquals (matrix.rows(), 0);
  validateEquals (matrix.cols(), 0);

  // Add columns
  matrix.append (DynamicVector<int>(3,1));
  validateEquals (matrix.size(), 3);
  validateEquals (matrix.rows(), 3);
  validateEquals (matrix.cols(), 1);
  validateEquals (matrix(0,0), 1);
  validateEquals (matrix(1,0), 1);
  validateEquals (matrix(2,0), 1);
  matrix.append (DynamicVector<int>(3,2));
  validateEquals (matrix.size(), 6);
  validateEquals (matrix.rows(), 3);
  validateEquals (matrix.cols(), 2);
  validateEquals (matrix(0,0), 1);
  validateEquals (matrix(1,0), 1);
  validateEquals (matrix(2,0), 1);
  validateEquals (matrix(0,1), 2);
  validateEquals (matrix(1,1), 2);
  validateEquals (matrix(2,1), 2);
  matrix.appendFront (DynamicVector<int>(3,3));
  validateEquals (matrix.size(), 9);
  validateEquals (matrix.rows(), 3);
  validateEquals (matrix.cols(), 3);
  validateEquals (matrix(0,0), 3);
  validateEquals (matrix(1,0), 3);
  validateEquals (matrix(2,0), 3);
  validateEquals (matrix(0,1), 1);
  validateEquals (matrix(1,1), 1);
  validateEquals (matrix(2,1), 1);
  validateEquals (matrix(0,2), 2);
  validateEquals (matrix(1,2), 2);
  validateEquals (matrix(2,2), 2);

  // Shift columns
  matrix.shiftSetFirst (DynamicVector<int>(3,4));
  validateEquals (matrix.size(), 9);
  validateEquals (matrix.rows(), 3);
  validateEquals (matrix.cols(), 3);
  validateEquals (matrix(0,0), 4);
  validateEquals (matrix(1,0), 4);
  validateEquals (matrix(2,0), 4);
  validateEquals (matrix(0,1), 3);
  validateEquals (matrix(1,1), 3);
  validateEquals (matrix(2,1), 3);
  validateEquals (matrix(0,2), 1);
  validateEquals (matrix(1,2), 1);
  validateEquals (matrix(2,2), 1);

  DynamicColumnMatrix<int> matrix2;
  matrix2.append (matrix);
  validateEquals (matrix2.size(), 9);
  validateEquals (matrix2.rows(), 3);
  validateEquals (matrix2.cols(), 3);
  validateEquals (matrix2(0,0), 4);
  validateEquals (matrix2(1,0), 4);
  validateEquals (matrix2(2,0), 4);
  validateEquals (matrix2(0,1), 3);
  validateEquals (matrix2(1,1), 3);
  validateEquals (matrix2(2,1), 3);
  validateEquals (matrix2(0,2), 1);
  validateEquals (matrix2(1,2), 1);
  validateEquals (matrix2(2,2), 1);

  // Remove columns from intermediate positions
  matrix2.remove (1);
  validateEquals (matrix2.size(), 6);
  validateEquals (matrix2.rows(), 3);
  validateEquals (matrix2.cols(), 2);
  validateEquals (matrix2(0,0), 4);
  validateEquals (matrix2(1,0), 4);
  validateEquals (matrix2(2,0), 4);
  validateEquals (matrix2(0,1), 1);
  validateEquals (matrix2(1,1), 1);
  validateEquals (matrix2(2,1), 1);
  matrix2.remove (0);
  validateEquals (matrix2.size(), 3);
  validateEquals (matrix2.rows(), 3);
  validateEquals (matrix2.cols(), 1);
  validateEquals (matrix2(0,0), 1);
  validateEquals (matrix2(1,0), 1);
  validateEquals (matrix2(2,0), 1);

  matrix2.clear ();
  validateEquals (matrix2.size(), 0);
  validateEquals (matrix2.rows(), 0);
  validateEquals (matrix2.cols(), 0);

}

}} // namespace tarch, la
