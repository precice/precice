#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "GramSchmidtTest.h"
#include "tarch/la/Matrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/GramSchmidt.h"
#include "tarch/la/MatrixMatrixOperations.h"
#include "tarch/la/TransposedMatrix.h"
#include "tarch/la/ScalarOperations.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(tarch::la::GramSchmidtTest)

namespace tarch {
namespace la {

GramSchmidtTest::GramSchmidtTest ()
:
  TestCase ("tarch::la::GramSchmidtTest")
{}

void GramSchmidtTest::run ()
{
  testMethod (testModifiedGramSchmidt);
}

void GramSchmidtTest::testModifiedGramSchmidt ()
{
  Matrix<4,4,double> A;
  DynamicMatrix<double> dynA (4, 4);

  // Set values according to Hilbert matrix.
  for (int i=0; i < 4; i++) {
     for (int j=0; j < 4; j++) {
        double entry = 1.0 / static_cast<double>(i + j + 1);
        A(i,j) = entry;
        dynA(i,j) = entry;
    }
  }

  Matrix<4,4,double> Acopy(A);
  Matrix<4,4,double> Q;
  Matrix<4,4,double> R;
  DynamicMatrix<double> dynAcopy (dynA);
  DynamicMatrix<double> dynQ (4, 4);
  DynamicMatrix<double> dynR (4, 4);
  modifiedGramSchmidt (Acopy, Q, R);
  modifiedGramSchmidt (dynAcopy, dynQ, dynR);

  Matrix<4,4,double> QTQ(0.0);
  DynamicMatrix<double> dynQTQ (4, 4, 0.0);
  multiply (transpose(Q), Q, QTQ);
  multiply (transpose(dynQ), dynQ, dynQTQ);

  for (int i=0; i < 4; i++) {
    for (int j=0; j < 4; j++) {
      if (i == j) {
        validate (equals(QTQ(i,j), 1.0, 1e-12));
        validate (equals(dynQTQ(i,j), 1.0, 1e-12));
      }
      else {
        validate (equals(QTQ(i,j), 0.0, 1e-12));
        validate (equals(dynQTQ(i,j), 0.0, 1e-12));
      }
    }
  }

  for (int i=0; i < 4; i++) {
     for (int j=0; j < 4; j++) {
        double value = 0.0;
        for ( int k=0; k < Q.cols(); k++ ) {
           value += Q(i,k) * R(k,j);
        }
        validate (equals(value, A(i,j)));
        value = 0.0;
     }
  }
}


}} // namespace tarch, la
