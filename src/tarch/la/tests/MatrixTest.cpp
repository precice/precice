#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#include "MatrixTest.h"
#include "tarch/la/Matrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicColumnMatrix.h"
#include "tarch/la/MatrixAssign.h"
#include "tarch/la/MatrixAssignList.h"
#include "tarch/la/MatrixOperations.h"
#include "tarch/la/MatrixMatrixOperations.h"
#include "tarch/la/TransposedMatrix.h"
#include <string>
#include <sstream>

#include "tarch/tests/TestCaseFactory.h"
registerTest(tarch::la::MatrixTest)

namespace tarch {
namespace la {

MatrixTest::MatrixTest ()
:
  TestCase ("tarch::la::MatrixTest")
{}

void MatrixTest::run ()
{
  testMethod (testConstruction);
  testMethod (testAssignment);
  testMethod (testMatrixOperations);
  testMethod (testMatrixMatrixOperations);
  testMethod (testTransposedMatrix);
}

void MatrixTest::testConstruction ()
{
  Matrix<1,2,int> matrix(1);
  validateEquals (matrix.size(), 2);
  validateEquals (matrix.rows(), 1);
  validateEquals (matrix.cols(), 2);
  validateEquals (matrix(0,0), 1);
  validateEquals (matrix(0,1), 1);

  Matrix<1,2,int> matrix2(matrix);
  validateEquals (matrix(0,0), 1);
  validateEquals (matrix(0,1), 1);

  DynamicMatrix<int> dynmatrix(1, 2, 1);
  validateEquals (dynmatrix.size(), 2);
  validateEquals (dynmatrix.rows(), 1);
  validateEquals (dynmatrix.cols(), 2);
  validateEquals (dynmatrix(0,0), 1);
  validateEquals (dynmatrix(0,1), 1);

  DynamicMatrix<int> dynmatrix2(dynmatrix);
  validateEquals (dynmatrix(0,0), 1);
  validateEquals (dynmatrix(0,1), 1);

  DynamicColumnMatrix<int> dyncolmatrix(1, 2, 1);
  validateEquals (dyncolmatrix.size(), 2);
  validateEquals (dyncolmatrix.rows(), 1);
  validateEquals (dyncolmatrix.cols(), 2);
  validateEquals (dyncolmatrix(0,0), 1);
  validateEquals (dyncolmatrix(0,1), 1);
}

void MatrixTest:: testAssignment ()
{
  Matrix<2,2,int> matrix(1);
  assignList(matrix) = 1, 2, 3, 4;
  validateEquals (matrix(0,0), 1);
  validateEquals (matrix(0,1), 2);
  validateEquals (matrix(1,0), 3);
  validateEquals (matrix(1,1), 4);

  DynamicMatrix<int> dynmatrix(2, 2, 2);
  validateEquals (dynmatrix(0,0), 2);
  validateEquals (dynmatrix(0,1), 2);
  validateEquals (dynmatrix(1,0), 2);
  validateEquals (dynmatrix(1,1), 2);

  assign(dynmatrix) = matrix;
  validateEquals (dynmatrix(0,0), 1);
  validateEquals (dynmatrix(0,1), 2);
  validateEquals (dynmatrix(1,0), 3);
  validateEquals (dynmatrix(1,1), 4);

  assign(dynmatrix) = 5;
  validateEquals (dynmatrix(0,0), 5);
  validateEquals (dynmatrix(0,1), 5);
  validateEquals (dynmatrix(1,0), 5);
  validateEquals (dynmatrix(1,1), 5);

  assign(matrix) = dynmatrix;
  validateEquals (matrix(0,0), 5);
  validateEquals (matrix(0,1), 5);
  validateEquals (matrix(1,0), 5);
  validateEquals (matrix(1,1), 5);

  Matrix<1,2,int> matrix2;
  assignList(matrix2) = 1, 2;
  validateEquals (matrix2(0,0), 1);
  validateEquals (matrix2(0,1), 2);
}

void MatrixTest::testMatrixOperations ()
{
  // Test computing determinant
  Matrix<3,3,int> matrix;
  assignList(matrix) = 1, 2, 3,
                       4, 5, 6,
                       7, 8, 9;
  validateEquals (0, det3x3(matrix));
  assignList(matrix) = -1, 2, 3,
                       4, -5, 2,
                       -2, 3, 1;
  validateEquals (1, det3x3(matrix)); // computed by octave

  // Test streaming
  Matrix<2,2,int> matrix2;
  assignList(matrix2) = 1, 2, 3, 4;
  std::ostringstream stream;
  stream << matrix2;
  validateEquals (stream.str(), std::string("1, 2; 3, 4"));
  // Test matrix multiply scalar
  matrix2 =matrix2*2;
  validateEquals(matrix2(0,0),2);
  validateEquals(matrix2(0,1),4);
  validateEquals(matrix2(1,0),6);
  validateEquals(matrix2(1,1),8);
  // Test matrix add matrix
//  Matrix<2,2,int> matrix3;
//  assignList(matrix3) = 1, 2, 3, 4;
//  matrix3=matrix3+matrix2;
//  validateEquals(matrix3(0,0),3);
//  validateEquals(matrix3(0,1),6);
//  validateEquals(matrix3(1,0),9);
//  validateEquals(matrix3(1,1),12);
  // Test matrix square
  Matrix<2,2,double> matrix4;
  assignList(matrix4) = 4.0, 9.0, 16.0, 25.0;
  matrix4=sqrt(matrix4);
  validateEquals(matrix4(0,0),2.0);
  validateEquals(matrix4(0,1),3.0);
  validateEquals(matrix4(1,0),4.0);
  validateEquals(matrix4(1,1),5.0);
}

void MatrixTest::testMatrixMatrixOperations ()
{
  Matrix<2,3,int> lMatrix;
  Matrix<3,2,int> rMatrix;
  Matrix<2,2,int> result(0);
  assignList(lMatrix) = 1, 2, 3,
                        4, 5, 6;
  assignList(rMatrix) = 6, 5,
                        4, 3,
                        2, 1;

  // Matrix matrix multiplication
  multiply (lMatrix, rMatrix, result);
  validateEquals (result(0,0), 20);
  validateEquals (result(0,1), 14);
  validateEquals (result(1,0), 56);
  validateEquals (result(1,1), 41);

  // Bitwise comparison
  Matrix<2,3,int> matrixA(1);
  Matrix<2,3,int> matrixB(2);
  validate (matrixA == matrixA);
  validate (! (matrixA == matrixB));

  // Test equalsReturnIndex
  Matrix<2,3,double> matrix1(1);
  Matrix<2,3,double> matrix2(2);
  int i=equalsReturnIndex(matrix1,matrix2);
  validateEquals(i,0);
}

void MatrixTest:: testTransposedMatrix ()
{
  Matrix<3,2,int> matrix;
  assignList(matrix) = 1, 2,
                       3, 4,
                       5, 6;
  typedef Matrix<3,2,int> Matrix;
  TransposedMatrix<Matrix>& transposed = transpose(matrix);
  validateEquals (transposed(0,0), 1);
  validateEquals (transposed(0,1), 3);
  validateEquals (transposed(0,2), 5);
  validateEquals (transposed(1,0), 2);
  validateEquals (transposed(1,1), 4);
  validateEquals (transposed(1,2), 6);

  validate (transpose(transposed) == matrix);
}

}} // namespace tarch, la
