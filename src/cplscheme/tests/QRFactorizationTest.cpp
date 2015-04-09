#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#include "QRFactorizationTest.h"
#include "cplscheme/impl/QRFactorization.hpp"
#include <Eigen/Dense>
#include "tarch/la/Matrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/GramSchmidt.h"
#include "tarch/la/MatrixMatrixOperations.h"
#include "tarch/la/TransposedMatrix.h"
#include "tarch/la/ScalarOperations.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::QRFactorizationTest)

namespace precice {
namespace cplscheme {
namespace tests {

QRFactorizationTest::QRFactorizationTest ()
:
  TestCase ("precice::cplscheme::QRFactorizationTest")
{}

void QRFactorizationTest::run ()
{
  testMethod (testQRFactorization);
}

void QRFactorizationTest::testQRFactorization ()
{
  int m = 6, n = 8;
  Eigen::MatrixXd A(n,m);
  tarch::la::DynamicMatrix<double> dynA (n, m);
  
  // Set values according to Hilbert matrix.
  for (int i=0; i < n; i++) {
     for (int j=0; j < m; j++) {
        double entry = 1.0 / static_cast<double>(i + j + 1);
        A(i,j) = entry;
        dynA(i,j) = entry;
    }
  }

  // compute QR factorization of A via successive inserting of columns
  impl::QRFactorization qr_1(A);
  
  tarch::la::DynamicMatrix<double> dynAcopy (dynA);
  tarch::la::DynamicMatrix<double> dynQ (n, m);
  tarch::la::DynamicMatrix<double> dynR (m, m);
  
  // compute QR factorization of A via modifiedGramSchmidt en block
  tarch::la::modifiedGramSchmidt (dynAcopy, dynQ, dynR);

  // test if Q^TQ equals identity 
  testQTQequalsIdentity(qr_1.matrixQ());
  testQTQequalsIdentity(dynQ);
  // test if QR equals A
  testQRequalsA(qr_1.matrixQ(), qr_1.matrixR(), A);
  
  /**
   * *************** deleting/adding Columns ************************
   */
  
  // ----------- delete last column ---------------
  Eigen::MatrixXd A_prime1 = A;
  Eigen::VectorXd col6 = A.col(m-1);
  A_prime1.conservativeResize(n,m-1);
  qr_1.deleteColumn(m-1);
  
  // A_prime1 = A(1:n, 1:m-1)
  
  // test if Q^TQ equals identity 
  testQTQequalsIdentity(qr_1.matrixQ());
  // test if QR equals A
  testQRequalsA(qr_1.matrixQ(), qr_1.matrixR(), A_prime1);
  
  // ----------- delete first column ---------------
  Eigen::MatrixXd A_prime2 = A_prime1;
  Eigen::VectorXd col1 = A.col(0);
  for(int i=0; i<A_prime2.rows(); i++)
    for(int j=0; j<A_prime2.cols()-1; j++)
      A_prime2(i,j) = A_prime2(i,j+1);
  A_prime2.conservativeResize(n,A_prime2.cols()-1);
  qr_1.deleteColumn(0);
  
  // A_prime2 = A(1:n, 2:m-1)
  
  // test if Q^TQ equals identity 
  testQTQequalsIdentity(qr_1.matrixQ());
  // test if QR equals A
  testQRequalsA(qr_1.matrixQ(), qr_1.matrixR(), A_prime2);
  
  // ----------- add first column -----------------
  qr_1.insertColumn(0, col1); 
  // test if Q^TQ equals identity 
  testQTQequalsIdentity(qr_1.matrixQ());
  // test if QR equals A
  testQRequalsA(qr_1.matrixQ(), qr_1.matrixR(), A_prime1);
  
  // ----------- add last column -----------------
  qr_1.insertColumn(qr_1.cols()-1, col6); // ?? 
  // test if Q^TQ equals identity 
  testQTQequalsIdentity(qr_1.matrixQ());
  // test if QR equals A
  testQRequalsA(qr_1.matrixQ(), qr_1.matrixR(), A);
  
  // ----------- delete middle column ---------------
  int k = 3;
  Eigen::MatrixXd A_prime3 = A;
  Eigen::VectorXd colk = A.col(k);
  for(int i=k; i<A_prime3.rows(); i++)
    for(int j=0; j<A_prime3.cols()-1; j++)
      A_prime3(i,j) = A_prime3(i,j+1);
  A_prime3.conservativeResize(n,A_prime3.cols()-1);
  qr_1.deleteColumn(k);
  // test if Q^TQ equals identity 
  testQTQequalsIdentity(qr_1.matrixQ());
  // test if QR equals A
  testQRequalsA(qr_1.matrixQ(), qr_1.matrixR(), A_prime3);
  
  // ----------- add middle column -----------------
  qr_1.insertColumn(k, colk);
  // test if Q^TQ equals identity 
  testQTQequalsIdentity(qr_1.matrixQ());
  // test if QR equals A
  testQRequalsA(qr_1.matrixQ(), qr_1.matrixR(), A);
  
}


void QRFactorizationTest::testQRequalsA(
  Eigen::MatrixXd& Q, 
  Eigen::MatrixXd& R, 
  Eigen::MatrixXd& A)
{
  Eigen::MatrixXd A_prime = Q*R;
  for (int i=0; i < A.rows(); i++) {
     for (int j=0; j < A.cols(); j++) {
        validate (tarch::la::equals(A_prime(i,j), A(i,j)));
     }
  }
}


void QRFactorizationTest::testQTQequalsIdentity(
  Eigen::MatrixXd& Q)
{
  Eigen::MatrixXd QTQ = Q.transpose() * Q;
  // test if Q^TQ equals identity 
  for (int i=0; i < QTQ.rows(); i++) {
    for (int j=0; j < QTQ.cols(); j++) {
      if (i == j) {
        validate (tarch::la::equals(QTQ(i,j), 1.0, 1e-12));
      }
      else {
        validate (tarch::la::equals(QTQ(i,j), 0.0, 1e-12));
      }
    }
  }
}

void QRFactorizationTest::testQTQequalsIdentity(
  tarch::la::DynamicMatrix< double >& dynQ)
{
  tarch::la::DynamicMatrix<double> dynQTQ (dynQ.cols(), dynQ.cols(), 0.0);
  // compute Q^TQ for modifiedGramSchmidt en block version
  tarch::la::multiply (tarch::la::transpose(dynQ), dynQ, dynQTQ);
  
  // test if Q^TQ equals identity 
  for (int i=0; i < dynQTQ.rows(); i++) {
    for (int j=0; j < dynQTQ.cols(); j++) {
      if (i == j) {
        validate (tarch::la::equals(dynQTQ(i,j), 1.0, 1e-12));
      }
      else {
        validate (tarch::la::equals(dynQTQ(i,j), 0.0, 1e-12));
      }
    }
  }
}




}}} // namespace precice, cplscheme
