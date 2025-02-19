#pragma once
/*
 * SVDFactorization.hpp
 *
 *  Created on: Feb 5, 2016
 *      Author: Klaudius Scheufele
 */

#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>

#include "acceleration/impl/ParallelMatrixOperations.hpp"
#include "acceleration/impl/Preconditioner.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"

// ------- CLASS DEFINITION

namespace precice {
namespace acceleration {
namespace impl {

/**
 * @brief Class that provides functionality to maintain a SVD decomposition of a matrix
 * via successive rank-1 updates and truncation with respect to the truncation threshold eps.
 */
class SVDFactorization {
public:
  // Eigen
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::VectorXd Vector;

  /**
   * @brief Constructor.
   */
  SVDFactorization(
      double            eps,
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
  template <typename Derived1, typename Derived2>
  void update(
      const Eigen::MatrixBase<Derived1> &A,
      const Eigen::MatrixBase<Derived2> &B)
  {
    PRECICE_TRACE();
    PRECICE_ASSERT(_initialized);
    /** updates the truncated svd factorization of the Jacobian with a rank-1 modification
     *
     * \psi * \sigma * \phi <-- \psi * \sigma * \phi + A * B^T
     *
     */
    if (_initialSVD) {
      PRECICE_ASSERT(A.rows() == _rowsA, A.rows(), _rowsA);
      PRECICE_ASSERT(B.rows() == _rowsB, B.rows(), _rowsB);
    } else {
      PRECICE_ASSERT(A.cols() == B.cols(), A.cols(), B.cols());
      _rowsA = A.rows();
      _rowsB = B.rows();
      _cols  = 0;
      _psi   = Matrix::Zero(_rowsA, 0);
      _phi   = Matrix::Zero(_rowsB, 0);
      _sigma = Vector::Zero(0);
    }

    /** (1): compute orthogonal basis P of (I-\psi\psi^T)A
     */
    Matrix Atil(_psi.cols(), A.cols()); // Atil is of size (K_bar x m)

    // Atil := \psi^T *A
    // local computation of \psi^T * A and allreduce_sum to Atil (global), stored local on each proc
    _parMatrixOps->multiply(_psi.transpose(), A, Atil, (int) _psi.cols(), _globalRowsA, (int) A.cols());

    // Ptil := (I-\psi\psi^T)A
    // Atil is local on each proc, thus fully local computation, embarrassingly parallel
    Matrix Ptil = A - _psi * Atil;

    // compute orthogonal basis P of Ptil, i.e., QR-dec (P, R_A) = QR(Ptil)
    Matrix P, R_A;
    computeQRdecomposition(Ptil, _globalRowsA, P, R_A);

    /**  (2): compute orthogonal basis Q of (I-\phi\phi^T)B
     */
    Matrix Btil(_phi.cols(), B.cols()); // Btil is of size (K_bar x m)
    // Btil := \phi^T *B
    _parMatrixOps->multiply(_phi.transpose(), B, Btil, (int) _phi.cols(), _globalRowsB, (int) B.cols());
    // Qtil := (I-\phi\phi^T)B
    Matrix Qtil = B - _phi * Btil;

    // compute orthogonal basis Q of Qtil, i.e., QR-dec (Q, R_B) = QR(Qtil)
    Matrix Q, R_B;
    computeQRdecomposition(Qtil, _globalRowsB, Q, R_B);

    //     e_orthModes.stop(true);

    /** (3) construct matrix \f$ K \in (K_bar + m -x) x (K_bar +m -y) \f$ if
     *      x .. deleted columns in P -> (m-x) new modes from A (rows of R_A)
     *      y .. deleted columns in Q -> (m-y) new modes from B (rows of R_B)
     *
     *      [ \sigma  0] + [ Atil ] * [ Btil ]^T
     *      [    0    0]   [ R_A  ]   [ R_B  ]
     *  (stored local on each proc).
     */
    Matrix K = Matrix::Zero(_psi.cols() + R_A.rows(), _psi.cols() + R_B.rows());
    Matrix K_A(_psi.cols() + R_A.rows(), Atil.cols());
    Matrix K_B(_phi.cols() + R_B.rows(), Btil.cols());

    for (int i = 0; i < _sigma.size(); i++)
      K(i, i) = _sigma(i);

    K_A.block(0, 0, Atil.rows(), Atil.cols())         = Atil;
    K_A.block(Atil.rows(), 0, R_A.rows(), R_A.cols()) = R_A;
    K_B.block(0, 0, Btil.rows(), Btil.cols())         = Btil;
    K_B.block(Btil.rows(), 0, R_B.rows(), R_B.cols()) = R_B;
    K += K_A * K_B.transpose();

    // compute svd of K
    Eigen::JacobiSVD<Matrix> svd(K, Eigen::ComputeThinU | Eigen::ComputeThinV);
    _sigma         = svd.singularValues();
    auto &psiPrime = svd.matrixU();
    auto &phiPrime = svd.matrixV();
    //     e_matK.stop(true);

    /** (4) rotate left and right subspaces
     */
    Matrix rotLeft(_rowsA, _psi.cols() + P.cols());
    Matrix rotRight(_rowsB, _phi.cols() + Q.cols());

    rotLeft.block(0, 0, _rowsA, _psi.cols())         = _psi;
    rotLeft.block(0, _psi.cols(), _rowsA, P.cols())  = P;
    rotRight.block(0, 0, _rowsB, _phi.cols())        = _phi;
    rotRight.block(0, _phi.cols(), _rowsB, Q.cols()) = Q;

    // [\psi,P] is distributed block-row wise, but \psiPrime is local on each proc, hence local mult.
    _psi = rotLeft * psiPrime;
    _phi = rotRight * phiPrime;

    //     e_rot.stop(true);

    /** (5) truncation of SVD
     */
    _cols = _sigma.size();

    int waste = 0;
    for (int i = 0; i < (int) _sigma.size(); i++) {
      if (_sigma(i) < _sigma(0) * _truncationEps) {
        _cols = i;
        waste = _sigma.size() - i;
        break;
      }
    }
    _waste += waste;

    _psi.conservativeResize(_rowsA, _cols);
    _phi.conservativeResize(_rowsB, _cols);
    _sigma.conservativeResize(_cols);
    PRECICE_ASSERT(_sigma(0) >= 0.0);
    PRECICE_DEBUG("SVD factorization of Jacobian is truncated to {} DOFs. Cut off {} DOFs", _cols, waste);

    _initialSVD = true;
  }

  /**
   * @brief: initializes the updated SVD factorization, i.e., sets the object for
   * parallel matrix-matrix operations and the number of global rows.
   */
  void initialize(PtrParMatrixOps parMatOps, int globalRowsA, int globalRowsB);

  /**
   * @brief: resets the SVD factorization
   */
  void reset();

  /**
   * @brief: returns a matrix representation of the orthogonal matrix Psi, A = Psi * Sigma * Phi^bs
   */
  Matrix &matrixPsi();

  /**
   * @brief: returns a matrix representation of the orthogonal matrix Sigma, A = Psi * Sigma * Phi^T
   */
  Vector &singularValues();

  /**
   * @brief: returns a matrix representation of the orthogonal matrix Phi, A = Psi * Sigma * Phi^T
   */
  Matrix &matrixPhi();

  /// @brief: returns the number of columns in the QR-decomposition
  int cols();

  /// @brief: returns the rank of the truncated SVD factorization
  Rank rank();

  /// @brief: returns the total number of truncated modes since last call to this method
  int getWaste();

  /// @brief: sets the threshold for the truncation of the SVD factorization
  void setThreshold(double eps);

  /// @brief: returns the truncation threshold for the SVD
  double getThreshold();

  /// @brief: applies the preconditioner to the factorized and truncated representation of the Jacobian matrix
  // void applyPreconditioner();

  /// @brief: appplies the inverse preconditioner to the factorized and truncated representation of the Jacobian matrix
  // void revertPreconditioner();

  void setPrecondApplied(bool b);

  /// @brief: enables or disables an additional QR-2 filter for the QR-decomposition
  void setApplyFilterQR(bool b, double eps = 1e-3);

  // bool isPrecondApplied();

  bool isSVDinitialized();

  /// Optional file-stream for logging output
  void setfstream(std::fstream *stream);

private:
  /** @brief: computes the QR decomposition of a matrix A of type A = PSI^T*A \in R^(rank x n)
   *
   *  This method computes a dedicated QR factorization [Q,R] = PSI^T*A.
   *  In case of linear dependence, i.e., if the next column cannot be orthogonalized to the
   *  columns in Q, this columns is deleted from Q, however not from R. To account for the
   *  missing column in Q, the corresponding line in R is deleted, though, the column in R
   *  is still filled with the projection weights, such that the dimension (cols) of R_A is
   *  aligned with the dimension (cols) of PSI^T*A, from which the matrix K is composed of.
   *
   *  The threshold parameter eps, indicates whether a column is seen to be in the column space
   *  of Q via the criterium ||v_orth|| / ||v|| <= eps (cmp. QR2 Filter)
   */
  void computeQRdecomposition(Matrix const &A, int globalRows, Matrix &Q, Matrix &R);

  logging::Logger _log{"acceleration::SVDFactorization"};

  /// @brief: preconditioner for least-squares system if vectorial system is used.
  PtrPreconditioner _preconditioner;

  /// @brief: object for parallel matrix operations, i.e., parallel mat-mat/ mat-vec multiplications
  PtrParMatrixOps _parMatrixOps = nullptr;

  /// @brief: SVD factorization of the matrix J = _psi * _sigma * _phi^T
  Matrix _psi;
  Matrix _phi;
  Vector _sigma;

  /// Number of rows of first multiplicator (on each proc, i.e., local)
  int _rowsA = 0;

  /// Number of rows of second multiplicator (on each proc, i.e., local)
  int _rowsB = 0;

  /// Number of columns, i.e., rank of the truncated svd
  int _cols = 0;

  /// Number of global rows of first multiplicator, i.e., sum of _rows for all procs
  int _globalRowsA = 0;

  /// Number of global rows of second multiplicator, i.e., sum of _rows for all procs
  int _globalRowsB = 0;

  // Total number of truncated modes after last call to method getWaste()
  int _waste = 0;

  /// Truncation parameter for the updated SVD decomposition
  double _truncationEps;

  /// Threshold for the QR2 filter for the QR decomposition.
  double _epsQR2 = 1e-3;

  /// true if the preconditioner has been applied appropriate to the updated SVD decomposition
  bool _preconditionerApplied = false;

  /// true, if ParallelMatrixOperations object is set, i.e., initialized
  bool _initialized = false;

  /// true, if at least one update has been made, i.e., the number of rows is known and a initial rank is given.
  bool _initialSVD = false;

  bool _applyFilterQR = false;
};

} // namespace impl
} // namespace acceleration
} // namespace precice

#endif /* PRECICE_NO_MPI */
