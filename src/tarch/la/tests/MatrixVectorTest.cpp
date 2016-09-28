#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "MatrixVectorTest.h"
#include "tarch/la/Matrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/MatrixAssignList.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/Vector.h"
#include "tarch/la/DynamicVector.h"
#include "tarch/la/VectorOperations.h"
#include "tarch/la/VectorAssign.h"
#include "tarch/la/VectorAssignList.h"
#include "tarch/la/MatrixVectorOperations.h"
#include <string>
#include <sstream>

#include "tarch/tests/TestCaseFactory.h"
registerTest(tarch::la::MatrixVectorTest)

namespace tarch {
namespace la {

MatrixVectorTest::MatrixVectorTest ()
:
  TestCase ("tarch::la::MatrixVectorTest")
{}

void MatrixVectorTest::run ()
{
  testMethod (testMultiplication);
  testMethod (testForwardSubstitution);
  testMethod (testBackSubstitution);
  testMethod (testSolveSystem3x3);
}

void MatrixVectorTest::testMultiplication ()
{
  Matrix<2,2,int> matrix;
  DynamicMatrix<int> dynmatrix(2,2);
  Vector<2,int> vector;
  DynamicVector<int> dynvector(2);
  assignList(matrix) = 1, 2, 3, 4;
  assignList(dynmatrix) = 1, 2, 3, 4;
  assignList(vector) = 1, 2;
  assignList(dynvector) = 1, 2;

  Vector<2,int> result(0);
  multiply (matrix, vector, result);
  validateEquals (result[0], 5);
  validateEquals (result[1], 11);

  assign(result) = 0;
  multiply (dynmatrix, vector, result);
  validateEquals (result[0], 5);
  validateEquals (result[1], 11);

  result = dynmatrix * vector;
  validateEquals (result[0], 5);
  validateEquals (result[1], 11);

  assign(result) = 0;
  multiply (dynmatrix, dynvector, result);
  validateEquals (result[0], 5);
  validateEquals (result[1], 11);

  result = dynmatrix * dynvector;
  validateEquals (result[0], 5);
  validateEquals (result[1], 11);

  assign(result) = 0;
  multiply (matrix, dynvector, result);
  validateEquals (result[0], 5);
  validateEquals (result[1], 11);

  result = matrix * dynvector;
  validateEquals (result[0], 5);
  validateEquals (result[1], 11);
}

void MatrixVectorTest::testForwardSubstitution ()
{
  DynamicMatrix<double> matrix (3, 3);
  assignList(matrix) =
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    2.0, 3.0, 0.0;
  DynamicVector<double> rhs(3);
  assignList(rhs) = 3.0, 2.0, 1.0;
  DynamicVector<double> x(3, 0.0);

  forwardSubstitution ( matrix, rhs, x );
  Vector<3,double> validX (3.0, -1.0, -2.0);
  validateWithParams1 ( equals(x, validX), x );
}

void MatrixVectorTest::testBackSubstitution ()
{
  DynamicMatrix<double> matrix (3, 3);
  assignList(matrix) = 1.0, 2.0, 3.0,
                       0.0, 4.0, 5.0,
                       0.0, 0.0, 6.0;
  DynamicVector<double> rhs (3);
  assignList(rhs) = 1.0, 2.0, 3.0;
  DynamicVector<double> x (3, 0.0);

  backSubstitution (matrix, rhs, x);

  validateNumericalEquals (x[2], 0.5 );
  validateNumericalEquals (x[1], -0.125 );
  validateNumericalEquals (x[0], -0.25 );
//
//  DynamicColumnMatrix<double> matrix2(5, 5);
//  assignList(matrix2) =  1.0, 2.0, 3.0, 4, 5,
//                         9.0, 8.0, 7.0, 6, 5,
//                         2.0, 4.0, 6.0, 8, 0,
//                         3.0, 4.0, 1.0, 2, 6,
//                         7.0, 6.0, 9.0, 8, 4;
//  DynamicVector<double> rhs2(5);
//  assignList(rhs2) = 1, 3, 5, 7, 8;
//  DynamicVector<double> x2(5, 0.0);
//  backSubstitution(matrix2, rhs, x);
//  validateNumericalEquals (x[4], 2);
}

void MatrixVectorTest:: testSolveSystem3x3 ()
{
  Matrix<3,3,double> matrix;
  assignList(matrix) = 1.0, 2.0, 3.0,
                       4.0, 5.0, 6.0,
                       2.0, 3.0, 2.0;
  Vector<3,double> rhs;
  assignList(rhs) = 1.0, 2.0, 3.0;
  Vector<3,double> x(0.0);

  rhs = solveSystem3x3 (matrix, rhs, x);
  validateNumericalEquals (x(0), -1.166666666666667);
  validateNumericalEquals (x(1),  2.333333333333333);
  validateNumericalEquals (x(2), -0.833333333333333);

  assignList(matrix) = 0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0,
                       1.0, 1.0, 1.0;
  assignList(rhs) = 0.5, 0.5, 1.0;
  assign(x) = 0.0;

  rhs = solveSystem3x3 (matrix, rhs, x);
  validateNumericalEquals (x(0),  0.0);
  validateNumericalEquals (x(1),  0.5);
  validateNumericalEquals (x(2),  0.5);
}


}} // namespace tarch, la
