#ifndef PRECICE_NO_MPI

#include "acceleration/impl/SVDFactorization.hpp"
#include <Eigen/Core>
#include <limits>
#include "utils/MasterSlave.hpp"

namespace precice {
namespace acceleration {
namespace impl {

SVDFactorization::SVDFactorization(
    double            eps,
    PtrPreconditioner preconditioner)
    : _preconditioner(preconditioner),
      _truncationEps(eps)
{
}

void SVDFactorization::initialize(
    PtrParMatrixOps parOps,
    int             globalRows)
{
  _parMatrixOps = parOps;
  _globalRows   = globalRows;
  _initialized  = true;
}

/*
void SVDFactorization::applyPreconditioner()
{
  PRECICE_TRACE();

  if(_psi.size() > 0 && _phi.size() > 0){
    // apply preconditioner: \psi_i * P_i, corresponds to Wtil_i * P_i, local!
    _preconditioner->apply(_psi);
    // apply preconditioner: \phi_i^T * P_i^{-1}, corresponds to Z_i * P_i^{-1}, local!
    // here, \phi^T should be preconditioned from right with inv_weights, i.e., the columns
    // of \phi^T are scaled. This is identical to scaling the rows of \phi, i.e., applying
    // P^{-1} * \phi
    _preconditioner->revert(_phi);
  }
  _preconditionerApplied = true;
}

void SVDFactorization::revertPreconditioner()
{
  PRECICE_TRACE();

  if(_psi.size() > 0 && _phi.size() > 0){
    // revert preconditioner: \psi_i * P_i^{-1}, corresponds to Wtil_i * P_i^{-1}, local!
    _preconditioner->revert(_psi);
    // revert preconditioner: \phi_i^T * P_i, corresponds to Z_i * P_i, local!
    // here, \phi^T should be preconditioned from right with _weights, i.e., the columns
    // of \phi^T are scaled. This is identical to scaling the rows of \phi, i.e., applying
    // P * \phi
    _preconditioner->apply(_phi);
  }
  _preconditionerApplied = false;
}
*/

void SVDFactorization::reset()
{
  _psi.resize(0, 0);
  _phi.resize(0, 0);
  _sigma.resize(0);
  _preconditionerApplied = false;
  _initialSVD            = false;
  _applyFilterQR         = false;
  _epsQR2                = 1e-3;
}

void SVDFactorization::computeQRdecomposition(
    Matrix const &A,
    Matrix &      Q,
    Matrix &      R)
{
  PRECICE_TRACE();

  // if nothing is linear dependent, the dimensions stay like this
  Q         = Matrix::Zero(A.rows(), A.cols());
  R         = Matrix::Zero(A.cols(), A.cols());
  int colsR = 0;
  int rowsR = 0;
  int colsQ = 0;

  // magic params:
  double omega = 0.;
  double theta = std::sqrt(2);

  // columns need to be inserted at the back, otherwise we would have to perform
  // givens rotations, to re-establish the upper diagonal form of R
  for (int colIndex = 0; colIndex < A.cols(); colIndex++) {

    // invariants:
    PRECICE_ASSERT(colsQ == rowsR, colsQ, rowsR);
    PRECICE_ASSERT(colsQ <= colIndex, colsQ, colIndex);

    Vector col = A.col(colIndex);

    // if system is quadratic; discard
    if (_globalRows == colIndex) {
      PRECICE_WARN("The matrix that is about to be factorized is quadratic, i.e., the new column cannot be orthogonalized; discard.");
      return;
    }

    /**
     * orthogonalize column "col" to columns in Q
     */
    Vector r        = Vector::Zero(R.rows()); // gram-schmidt coefficients orthogonalization
    Vector s        = Vector::Zero(R.rows()); // gram-schmidt coefficients re-orthogonalization
    Vector u        = Vector::Zero(A.rows()); // sum of projections
    double rho_orth = 0.;
    double rho0     = utils::MasterSlave::l2norm(col); // distributed l2norm;
    double rho00    = rho0;                            // save norm of col for QR2 filter crit.

    int  its            = 0;
    bool termination    = false;
    bool orthogonalized = true;
    // while col is not sufficiently orthogonal to Q
    while (!termination) {

      // take a gram-schmidt iteration
      u = Vector::Zero(A.rows());
      for (int j = 0; j < colsQ; j++) {

        // dot-product <_Q(:,j), v >
        Vector Qc = Q.col(j);
        // dot product <_Q(:,j), v> =: r_ij
        double r_ij = utils::MasterSlave::dot(Qc, col);
        // save r_ij in s(j) = column of R
        s(j) = r_ij;
        // u is the sum of projections r_ij * _Q(:,j) =  _Q(:,j) * <_Q(:,j), v>
        u += Q.col(j) * r_ij;
      }
      // add the gram-schmidt coefficients over all iterations of possible re-orthogonalizations
      r += s;
      // subtract projections from v, v is now orthogonal to columns of Q
      col -= u;

      // rho1 = norm of orthogonalized new column v_tilde (though not normalized)
      rho_orth = utils::MasterSlave::l2norm(col);
      // t = norm of _r(:,j) with j = colNum-1
      double norm_coefficients = utils::MasterSlave::l2norm(s); // distributed l2norm

      its++;

      // if ||v_orth|| is nearly zero, col is not well orthogonalized; discard
      if (rho_orth <= std::numeric_limits<double>::min()) {
        PRECICE_DEBUG("The norm of v_orthogonal is almost zero, i.e., failed to orthogonalize column v; discard.");
        orthogonalized = false;
        termination    = true;
      }

      /**   - test if reorthogonalization is necessary -
       *  rho0 = |v_init|, t = |r_(i,cols-1)|, rho_orth = |v_orth|
       *  rho_orth is small, if the new information incorporated in v is small,
       *  i.e., the part of v orthogonal to _Q is small.
       *  if rho_orth is very small it is possible, that we are adding (more or less)
       *  only round-off errors to the decomposition. Later normalization will scale
       *  this new information so that it is equally weighted as the columns in Q.
       *  To keep a good orthogonality, the gram-schmidt process is iterated for the
       *  bad orthogonalized column v_orth = col'
       *
       *  re-orthogonalize if: ||v_orth|| / ||v|| <= 1/theta
       */
      if (rho_orth * theta <= rho0 + omega * norm_coefficients) {
        // exit to fail if too many iterations
        if (its >= 4) {
          PRECICE_WARN("Matrix Q is not sufficiently orthogonal. Failed to rorthogonalize new column after 4 iterations. New column will be discarded.");
          orthogonalized = false;
          termination    = true;
        }

        // for re-orthogonalization
        rho0 = rho_orth;

      } else {
        termination = true;
      }
    }

    // if the QR2-filter crit. kicks in with threshold eps.
    if (_applyFilterQR && orthogonalized && rho_orth <= _epsQR2 * rho00) {
      orthogonalized = false;
    }

    // normalize col
    double rho = orthogonalized ? rho_orth : 1.0;
    col /= rho;
    r(rowsR) = rho;

    // as we always insert at the rightmost position, no need to shift
    // entries of R or apply givens rotations on the QR-dec to maintain
    // the upper triangular structure of R
    Q.col(colsQ) = col; // insert orthogonalized column to the right in Q
    R.col(colsR) = r;   // insert gram-schmidt coefficients to the right in R

    colsR++;
    PRECICE_ASSERT(colsR <= R.cols(), colsR, R.cols());
    rowsR++;
    PRECICE_ASSERT(rowsR <= R.rows(), rowsR, R.rows());
    colsQ++;

    // failed to orthogonalize the column, i.e., it is linear dependent;
    // modify the QR-dec such that it stays valid (column deleted) while
    // also staying aligned with the dimension of A, MUST have the same
    // number of cols (cannot delete from A)
    if (not orthogonalized) {

      colsQ--;
      PRECICE_ASSERT(colsQ >= 0, colsQ);
      rowsR--;
      PRECICE_ASSERT(rowsR >= 0, rowsR);
      // delete column that was just inserted (as it is not orthogonal to Q)
      Q.col(colsQ) = Vector::Zero(A.rows());
      // delete line in R that corresponds to the just inserted but not orthogonal column
      // as we always insert to the right, no shifting/ application of givens roatations is
      // necessary.
      // Note: The corresponding column from R with index colIndex is not deleted: dimensions must align with A.
      PRECICE_ASSERT(R(rowsR, colsR - 1) == 1.0, R(rowsR, colsR - 1));
      R.row(rowsR) = Vector::Zero(A.cols());
    }
  }
  // shrink matrices Q, R to actual size
  Q.conservativeResize(A.rows(), colsQ);
  R.conservativeResize(rowsR, colsR);
}

SVDFactorization::Matrix &SVDFactorization::matrixPhi()
{
  return _phi;
}

SVDFactorization::Matrix &SVDFactorization::matrixPsi()
{
  return _psi;
}

SVDFactorization::Vector &SVDFactorization::singularValues()
{
  return _sigma;
}

void SVDFactorization::setPrecondApplied(bool b)
{
  _preconditionerApplied = b;
}

void SVDFactorization::setApplyFilterQR(bool b, double eps)
{
  _applyFilterQR = b;
  _epsQR2        = eps;
}

/*
bool SVDFactorization::isPrecondApplied()
{
  return _preconditionerApplied;
}
*/

bool SVDFactorization::isSVDinitialized()
{
  return _initialSVD;
}

void SVDFactorization::setThreshold(double eps)
{
  _truncationEps = eps;
}

double SVDFactorization::getThreshold()
{
  return _truncationEps;
}

int SVDFactorization::getWaste()
{
  int r  = _waste;
  _waste = 0;
  return r;
}

int SVDFactorization::cols()
{
  return _cols;
}

int SVDFactorization::rows()
{
  return _rows;
}

int SVDFactorization::rank()
{
  return _cols;
}

} // namespace impl
} // namespace acceleration
} // namespace precice

#endif // PRECICE_NO_MPI
