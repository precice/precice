#include "mapping/RadialBasisFctSolver.hpp"

#include <Eigen/QR>
#include "precice/types.hpp"
#include "utils/Event.hpp"

namespace precice {
namespace mapping {

Eigen::VectorXd RadialBasisFctSolver::solveConservative(const Eigen::VectorXd &inputData, Polynomial polynomial) const
{
  // TODO: Avoid temporary allocations
  // Au is equal to the eta in our PETSc implementation
  PRECICE_ASSERT(inputData.size() == _matrixA.rows());
  Eigen::VectorXd Au = _matrixA.transpose() * inputData;
  PRECICE_ASSERT(Au.size() == _matrixA.cols());

  // mu in the PETSc implementation
  Eigen::VectorXd out = _qrMatrixC.solve(Au);

  if (polynomial == Polynomial::SEPARATE) {
    Eigen::VectorXd epsilon = _matrixV.transpose() * inputData;
    PRECICE_ASSERT(epsilon.size() == _matrixV.cols());

    // epsilon = Q^T * mu - epsilon (tau in the PETSc impl)
    epsilon -= _matrixQ.transpose() * out;
    PRECICE_ASSERT(epsilon.size() == _matrixQ.cols());

    // out  = out - solveTranspose tau (sigma in the PETSc impl)
#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
    out -= static_cast<Eigen::VectorXd>(_qrMatrixQ.transpose().solve(-epsilon));
#else
    // Backwards compatible version
    Eigen::VectorXd    sigma(_matrixQ.rows());
    const Eigen::Index nonzero_pivots = _qrMatrixQ.nonzeroPivots();

    if (nonzero_pivots == 0) {
      sigma.setZero();
    } else {
      Eigen::VectorXd c(_qrMatrixQ.colsPermutation().transpose() * (-epsilon));

      _qrMatrixQ.matrixQR().topLeftCorner(nonzero_pivots, nonzero_pivots).template triangularView<Eigen::Upper>().transpose().conjugate().solveInPlace(c.topRows(nonzero_pivots));

      sigma.topRows(nonzero_pivots) = c.topRows(nonzero_pivots);
      sigma.bottomRows(_qrMatrixQ.rows() - nonzero_pivots).setZero();

      sigma.applyOnTheLeft(_qrMatrixQ.householderQ().setLength(nonzero_pivots));
      out -= sigma;
    }
#endif
  }
  return out;
}

Eigen::VectorXd RadialBasisFctSolver::solveConsistent(Eigen::VectorXd &inputData, Polynomial polynomial) const
{
  Eigen::VectorXd res;
  // Solve polynomial QR and substract it form the input data
  if (polynomial == Polynomial::SEPARATE) {
    res = _qrMatrixQ.solve(inputData);
    inputData -= (_matrixQ * res);
  }

  // Integrated polynomial (and separated)
  PRECICE_ASSERT(inputData.size() == _matrixA.cols());
  Eigen::VectorXd p = _qrMatrixC.solve(inputData);
  PRECICE_ASSERT(p.size() == _matrixA.cols());
  Eigen::VectorXd out = _matrixA * p;

  // Add the polynomial part again for separated polynomial
  if (polynomial == Polynomial::SEPARATE) {
    out += (_matrixV * res);
  }
  return out;
}

void RadialBasisFctSolver::clear()
{
  _matrixA   = Eigen::MatrixXd();
  _qrMatrixC = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>();
}

const Eigen::MatrixXd &RadialBasisFctSolver::getEvaluationMatrix() const
{
  return _matrixA;
}
} // namespace mapping
} // namespace precice
