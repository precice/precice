/*
 * SVDFactorization.hpp
 *
 *  Created on: Feb 5, 2016
 *      Author: Klaudius Scheufele
 */

#ifndef PRECICE_NO_MPI

#ifndef SVDFACTORIZATION_HPP_
#define SVDFACTORIZATION_HPP_

#include "mesh/SharedPointer.hpp"
#include "SharedPointer.hpp"
#include "ParallelMatrixOperations.hpp"
#include "QRFactorization.hpp"
#include "Preconditioner.hpp"
#include "tarch/logging/Log.h"
#include "utils/MasterSlave.hpp"
#include "utils/EventTimings.hpp"
#include <Eigen/Dense>
#include <limits>
#include <fstream>



// ------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Class that provides functionality to maintain a SVD decomposition of a matrix
 * via succesive rank-1 updates and truncation with respect to the truncation threshold eps.
 */
class SVDFactorization
{
public:
  // Eigen
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::VectorXd Vector;



  /**
   * @brief Constructor.
   */
  SVDFactorization(
      double eps,
      PtrPreconditioner preconditioner);

   /**
    * @brief Destructor, empty.
    */
   virtual ~SVDFactorization() {}

   /** @brief: updates the SVD decomposition with the rank-1 update A*B^T, i.e.,
    *               _psi * _sigma * _phi^T + A*B^T
    *  and overrides the internal SVD representation. After the update, the SVD is
    *  truncated according to the threshold _truncationEps
    */
   template<typename Derived1, typename Derived2>
   void update(
       const Eigen::MatrixBase<Derived1>& A,
       const Eigen::MatrixBase<Derived2>& B)
   {
     preciceTrace(__func__);
     utils::Event e("SVD update", true, true);
     assertion(_initialized);
     assertion(_preconditionerApplied);
     /** updates the truncated svd factorization of the Jacobian with a rank-1 modification
      *
      * \psi * \sigma * \phi <-- \psi * \sigma * \phi + A * B^T
      *
      * PRECONDITION: The internal representation of the SVD is preconditioned correctly and also the
      * incoming matrices A and B.
      */
     if(_initialSVD){
       assertion2(A.rows() == _rows, A.rows(), _rows);
       assertion2(B.rows() == _rows, B.rows(), _rows);
     }else{
       assertion2(A.rows() == B.rows(), A.rows(), B.rows()); assertion2(A.cols() == B.cols(), A.cols(), B.cols());
       _rows = A.rows();
       _cols = 0;
       _psi = Matrix::Zero(_rows, 0);
       _phi = Matrix::Zero(_rows, 0);
       _sigma = Vector::Zero(0);
     }

     /** (1): compute orthogonal basis P of (I-\psi\psi^T)A
      */
     Matrix Atil(_psi.cols(), A.cols());    // Atil is of size (K_bar x m)

     // Atil := \psi^T *A
     // local computation of \psi^T * A and allreduce_sum to Atil (global), stored local on each proc
     _parMatrixOps->multiply(_psi.transpose(), A, Atil, (int)_psi.cols(), _globalRows, (int)A.cols());

     // Ptil := (I-\psi\psi^T)A
     // Atil is local on each proc, thus fully local computation, embarrassingly parallel
     Matrix Ptil = A - _psi * Atil;

     // compute orthogonal basis P of Ptil, i.e., QR-dec (P, R_A) = QR(Ptil)
     QRFactorization qrA(0);        // TODO: maybe add filter, currently hard coded to NO_FILTER
     qrA.reset(Ptil, _globalRows);  // builds QR-factorization of matrix Ptil
     auto& P   = qrA.matrixQ();
     auto& R_A = qrA.matrixR();

     /**  (2): compute orthogonal basis Q of (I-\phi\phi^T)B
      */
     Matrix Btil(_phi.cols(), B.cols());    // Btil is of size (K_bar x m)
     // Btil := \phi^T *B
     _parMatrixOps->multiply(_phi.transpose(), B, Btil, (int)_phi.cols(), _globalRows, (int)B.cols());
     // Qtil := (I-\phi\phi^T)B
     Matrix Qtil = B - _phi * Btil;

     // compute orthogonal basis Q of Qtil, i.e., QR-dec (Q, R_B) = QR(Qtil)
     QRFactorization qrB(0);        // TODO: maybe add filter, currently hard coded to NO_FILTER
     qrB.reset(Qtil, _globalRows);  // builds QR-factorization of matrix Ptil
     auto& Q   = qrB.matrixQ();
     auto& R_B = qrB.matrixR();

     /** (3) construct matrix K \in (K_bar + m) x (K_bar +m)
      *      [ \sigma  0] + [ Atil ] * [ Btil ]^T
      *      [    0    0]   [ R_A  ]   [ R_B  ]
      *  (stored local on each proc).
      */
     Matrix K = Matrix::Zero(_psi.cols() + Atil.cols(), _psi.cols() + Atil.cols());
     Matrix K_A(_psi.cols() + Atil.cols(), Atil.cols());
     Matrix K_B(_phi.cols() + Btil.cols(), Btil.cols());
     assertion2(Atil.cols() == Btil.cols(), Atil.cols(), Btil.cols()); //TODO: check if this must be the case, how do we multiply if not?
     for(int i = 0; i < _sigma.size(); i++)
       K(i,i) = _sigma(i);

     K_A.block(0,0,Atil.rows(),Atil.cols()) = Atil;
     K_A.block(Atil.rows(), 0, R_A.rows(), R_A.cols()) = R_A;
     K_B.block(0,0,Btil.rows(),Btil.cols()) = Btil;
     K_B.block(Btil.rows(), 0, R_B.rows(), R_B.cols()) = R_B;
     K += K_A * K_B.transpose();

     // compute svd of K
     Eigen::JacobiSVD<Matrix> svd(K, Eigen::ComputeThinU | Eigen::ComputeThinV);
     _sigma = svd.singularValues();
     auto& psiPrime = svd.matrixU();
     auto& phiPrime = svd.matrixV();

     /** (4) rotate left and right subspaces
      */
     Matrix rotLeft(_rows, _psi.cols() + P.cols());
     Matrix rotRight(_rows, _phi.cols() + Q.cols());

     rotLeft.block(0,0,_rows, _psi.cols()) = _psi;
     rotLeft.block(0,_psi.cols(),_rows, P.cols()) = P;
     rotRight.block(0,0,_rows, _phi.cols()) = _phi;
     rotRight.block(0,_phi.cols(),_rows, Q.cols()) = Q;

     // [\psi,P] is distributed block-row wise, but \psiPrime is local on each proc, hence local mult.
     _psi = rotLeft * psiPrime;
     _phi = rotRight * phiPrime;

     /** (5) truncation of SVD
      */
     _cols = _sigma.size();

     for(int i = 0; i < (int)_sigma.size(); i++){
       if(_sigma(i) < (int)_sigma(0) * _truncationEps){
         _cols = i;
         break;
       }
     }
     int waste = _sigma.size()-i;

     _psi.conservativeResize(_rows, _cols);
     _phi.conservativeResize(_rows, _cols);
     _sigma.conservativeResize(_cols);
     preciceDebug("SVD factorization of Jacobian is truncated to "<<_cols<<" DOFs. Cut off "<<waste<<" DOFs");

     _initialSVD = true;
   }


   /**
    * @brief: initializes the updated SVD factorization, i.e., sets the object for
    * parallel matrix-matrix operations and the number of global rows.
    */
   void initialize(PtrParMatrixOps parMatOps, int globalRows);

   /**
    * @brief: resets the SVD factorization
    */
   void reset();

   /**
    * @brief: returns a matrix representation of the orthogonal matrix Psi, A = Psi * Sigma * Phi^T
    */
   Matrix& matrixPsi();

   /**
    * @brief: returns a matrix representation of the orthogonal matrix Sigma, A = Psi * Sigma * Phi^T
    */
   Vector& singularValues();

   /**
    * @brief: returns a matrix representation of the orthogonal matrix Phi, A = Psi * Sigma * Phi^T
    */
   Matrix& matrixPhi();

   /// @brief: returns the number of columns in the QR-decomposition
   int cols();

   /// @brief: returns the number of rows in the QR-decomposition
   int rows();

   /// @brief: returns the rank of the truncated SVD factorization
   int rank();

   /// @brief: sets the threshold for the truncation of the SVD factorization
   void setThreshold(double eps);

   /// @brief: applies the preconditioner to the factorized and truncated representation of the Jacobian matrix
   void applyPreconditioner();

   /// @brief: appplies the inverse preconditioner to the factorized and truncated representation of the Jacobian matrix
   void revertPreconditioner();

   void setPrecondApplied(bool b);

   bool isPrecondApplied();

   bool isSVDinitialized();

   // @brief optional file-stream for logging output
   void setfstream(std::fstream* stream);

private:


  /// @brief: Logging device.
  static tarch::logging::Log _log;

  /// @brief: preconditioner for least-squares system if vectorial system is used.
  PtrPreconditioner _preconditioner;

  /// @brief: object for parallel matrix operations, i.e., parallel mat-mat/ mat-vec multiplications
  PtrParMatrixOps _parMatrixOps;

  /// @brief: SVD factorization of the matrix J = _psi * _sigma * _phi^T
  Matrix _psi;
  Matrix _phi;
  Vector _sigma;

  /// @brief: number of rows (on each proc, i.e., local)
  int _rows;

  /// @brief number of columns, i.e., rank of the truncated svd
  int _cols;

  /// @brief: numer of global rows, i.e., sum of _rows for all procs
  int _globalRows;

  ///@brief: Truncation parameter for the updated SVD decomposition
  double _truncationEps;

  /// @brief: true if the preconditioner has been applied appropriate to the updated SVD decomposition
  bool   _preconditionerApplied;

  /// @brief: true, if ParallelMatrixOperations object is set, i.e., initialized
  bool _initialized;

  /// @brief: true, if at least one update has been made, i.e., the number of rows is known and a initial rank is given.
  bool _initialSVD;

  // @brief optional infostream that writes information to file
  std::fstream* _infostream;
  bool _fstream_set;

};

}}} // namespace precice, cplscheme, impl



#endif /* SVDFACTORIZATION_HPP_ */
#endif /* PRECICE_NO_MPI */
