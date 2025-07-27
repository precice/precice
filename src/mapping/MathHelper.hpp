#pragma once

#include <Eigen/Core>

namespace precice::utils {
/**
 * @brief Implements an iterative block scheme to determine the inverse of a lower triangular matrix
 * which is more efficient than Eigen's own trsm
 *
 * @param L the lower triangular matrix
 * @return Eigen::MatrixXd the inverse of L
 */
template <typename T>
inline Eigen::Matrix<T, -1, -1> invertLowerTriangularBlockwise(const Eigen::Matrix<T, -1, -1> &L)
{
  const int n = L.rows();

  // We copy over Eigen's block size heuristic to set the block size
  int blockSize = static_cast<int>(n / 8);
  blockSize     = (blockSize / 16) * 16;
  blockSize     = std::min(std::max(blockSize, 8), 128);

  // Handle small matrices by setting blockSize to n if it's zero
  if (blockSize == 0) {
    blockSize = n;
  }
  const int numBlocks = static_cast<int>((n + blockSize - 1) / blockSize);

  // Initialize the inverse matrix as a identity (enables inPlaceSolve)
  // TODO: maybe consider to return by reference, as Eigen doesn't recommend to return by value
  Eigen::Matrix<T, -1, -1> L_inv = Eigen::Matrix<T, -1, -1>::Identity(n, n);

  // Now we iterate over all blocks...
  for (int i = 0; i < numBlocks; ++i) {
    int i_start     = i * blockSize;
    int i_end       = std::min(i_start + blockSize, n);
    int i_blockSize = i_end - i_start;

    // Step 1: Invert the diagonal block
    L.block(i_start, i_start, i_blockSize, i_blockSize)
        .template triangularView<Eigen::Lower>()
        .solveInPlace(L_inv.block(i_start, i_start, i_blockSize, i_blockSize));

    // Step 2: Compute off-diagonal blocks
    for (int j = 0; j < i; ++j) {
      int j_start     = j * blockSize;
      int j_end       = std::min(j_start + blockSize, n);
      int j_blockSize = j_end - j_start;

      // computation for the off-diagonal block using the block Lij and Ljj (from L_inv)
      Eigen::Matrix<T, -1, -1> offDiagBlock = L.block(i_start, j_start, i_blockSize, j_blockSize) * L_inv.block(j_start, j_start, j_blockSize, j_blockSize)
                                                                                                        .template triangularView<Eigen::Lower>();
      // Accumulate the sum of L_ik * L_inv(k, j) for k = j+1 to i-1
      for (int k = j + 1; k < i; ++k) {
        int k_start     = k * blockSize;
        int k_end       = std::min(k_start + blockSize, n);
        int k_blockSize = k_end - k_start;

        // Update offDiagBlock with L_ik + L_inv_kj
        offDiagBlock += L.block(i_start, k_start, i_blockSize, k_blockSize) * L_inv.block(k_start, j_start, k_blockSize, j_blockSize);
      }

      // Fill in the off-diagonal block into the matrix
      // Note that, although L_inv is on both sides, there source and destination are different in memory
      L_inv.block(i_start, j_start, i_blockSize, j_blockSize) = -L_inv.block(i_start, i_start, i_blockSize, i_blockSize) * offDiagBlock;
    }
  }
  return L_inv;
}

/**
 * @brief For C = LL^T, compute the diagonal entries of the inverse kernel matrix.
 * @param decMatrixC Cholesky decomposition of the kernel matrix.
 */
inline Eigen::VectorXd computeInverseDiagonal(const Eigen::LLT<Eigen::MatrixXd> &decMatrixC)
{
  // 1. Compute the diagonal entries of the inverse kernel matrix:
  // We already have the Cholesky decomposition. So instead of solving for the
  // kernel matrix directly, we invert the lower triangular matrix of the
  // decomposition:
  // using A^{-1} = (L^T)^{-1}L^{-1} enables the computation of the diagonal
  // entries of A^{-1} by evaluating the product above:
  // A^{-1}_{ii} = sum_{k=1}^n L^{-T}_{ik} L^{-1}_{ki}

  // 1a: Compute the inverse of the lower triangular matrix L
  // Eigen::MatrixXd L_inv = L.inverse(); is not supported by Eigen (linker errors)
  // However, Eigen provides triangular solver (LAPACK::trsm), which can be used
  // to solve L * Linv = I
  // Eigen::MatrixXd L_inv = Eigen::MatrixXd::Identity(inSize, inSize);
  // _decMatrixC.matrixL().solveInPlace(L_inv);
  // which yields cubic complexity (BLAS level 3).

  // Here, we use our own implementation to compute the triangular inverse
  // (similar to LAPACK:trtri) more efficient, given that the RHS is also
  // triangular was unfortunately slower.

  // Solve L * Linv = I
  // Eigen::MatrixXd L_inv = Eigen::MatrixXd::Identity(n, n);
  // decMatrixC.matrixL().solveInPlace(L_inv);
  Eigen::MatrixXd L_inv = utils::invertLowerTriangularBlockwise<double>(decMatrixC.matrixL());

  // 1b: Compute the diagonal elements of A^{-1} by evaluating (L^T)^{-1}L^{-1}
  Eigen::VectorXd inverseDiagonal = (L_inv.array().square().colwise().sum()).transpose();

  return inverseDiagonal;
}

/**
 * @brief For C = QR, compute the diagonal entries of the inverse kernel matrix.
 * @param decMatrixC QR decomposition of the kernel matrix.
 */
inline Eigen::VectorXd computeInverseDiagonal(const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> &decMatrixC)
{
  // 1. Compute the diagonal entries of the inverse kernel matrix:
  // We could use
  //     diag_inv_A= _decMatrixC.inverse().diagonal();
  // as Eigen offers this for the QR decomposition, but the inverse() call is
  // more expensive than it has to be (as we compute all entries of the inverse). On the
  // other hand, using non strictily-positive definite functions is less relevant anyway.
  // Still, let's try the following:

  // We already have the QR decomposition. So instead of solving for the
  // kernel matrix directly, make use of the following:
  // A^{-1} = R^{-1}Q^{-1}
  // Since Q is orthogonal, Q^{-1} = Q^T
  // R is upper triangular and we need to compute the inverse (using backwards substitution)

  // enables the computation of the diagonal
  // entries of A^{-1} by evaluating the product above:
  // A^{-1} = R^{-1} Q^T
  // A^{-1}_{ii} = sum_{k=1}^n R^{-1}_{ik} Q^{T}_{ki}

  Eigen::VectorXd    inverseDiagonal;
  const Eigen::Index n = decMatrixC.matrixR().cols();

  // 1a: Compute the inverse of the lower triangular matrix L
  // Solve R * Rinv = I
  Eigen::MatrixXd R_inv = Eigen::MatrixXd::Identity(n, n);
  decMatrixC.matrixR().template triangularView<Eigen::Upper>().solveInPlace(R_inv);
  Eigen::VectorXi P = decMatrixC.colsPermutation().indices();
  Eigen::MatrixXd Q = decMatrixC.householderQ();

  // Now evaluate, caution with the column permutation
  inverseDiagonal.resize(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    inverseDiagonal(P(i)) = R_inv.row(i) * Q.transpose().col(P(i));
  }
  return inverseDiagonal;
}

/**
 * @brief Computes the Leave-One-Out Cross Validation error from a Cholesky decomposition.
 * @param choleskyDec Cholesky decomposition.
 * @param inputData Right hand side data of the linear system.
 * @return LOOCV error for a valid Cholesky decomposition and NaN otherwise.
 *
 * Implementation of LOOCV according to Rippa(1999), DOI: 10.1023/a:1018975909870
 */
template <typename DecompositionType>
double computeRippaLOOCVerror(const DecompositionType &choleskyDec, const Eigen::VectorXd &inputData)
{
  static_assert(std::is_same_v<DecompositionType, Eigen::LLT<Eigen::MatrixXd>> || std::is_same_v<DecompositionType, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>,
                "computeRippaLOOCVerror() only allows DecompositionType to be Eigen::LLT<Eigen::MatrixXd>> or Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>");

  if (choleskyDec.info() != Eigen::ComputationInfo::Success) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const Eigen::Index    n      = inputData.size();
  const Eigen::VectorXd lambda = choleskyDec.solve(inputData);

  const double loocv = std::sqrt((lambda.array() / computeInverseDiagonal(choleskyDec).array()).array().square().sum() / n);

  return loocv;
}

/**
 * @brief Computes an approximation such that 1/cond(LL^T) > returned result.
 *
 * Implementation based on the diagonal entries of the Cholesky decomposition matrix.
 * See also: "A Survey of Condition Number Estimation for Triangular Matrices" DOI: 10.1137/1029112
 *
 * @return approx < 1/cond(LL^T) and 0 if the decomposition is invalid.
 */
inline double approximateReciprocalConditionNumber(const Eigen::LLT<Eigen::MatrixXd> &choleskyDec)
{
  if (choleskyDec.info() != Eigen::ComputationInfo::Success) {
    return 0;
  }
  const Eigen::Index n = choleskyDec.matrixL().rows();

  double max_l = std::numeric_limits<double>::min();
  double min_l = std::numeric_limits<double>::max();

  for (Eigen::Index i = 0; i < n; i++) {
    const double lii = choleskyDec.matrixL()(i, i);
    min_l            = std::min(lii, min_l);
    max_l            = std::max(lii, max_l);
  }
  double rcond = min_l * min_l / (max_l * max_l);
  if (rcond < 0 || std::isnan(rcond))
    rcond = 0;

  return rcond;
}

/**
 * @brief Computes an approximation such that 1/cond(QR) >= returned result.
 *
 * Implementation based on the diagonal entries of the triangular R matrix of the QR decomposition.
 * See also: "A Survey of Condition Number Estimation for Triangular Matrices" https://doi.org/10.1137/1029112
 *
 * @return >= 1/cond(QR)
 */
inline double approximateReciprocalConditionNumber(const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> &qrDec)
{
  if (qrDec.info() != Eigen::ComputationInfo::Success) {
    return 0;
  }
  const Eigen::Index n = qrDec.matrixR().rows();

  double max_r = std::numeric_limits<double>::min();
  double min_r = std::numeric_limits<double>::max();

  for (Eigen::Index i = 0; i < n; i++) {
    const double rii = qrDec.matrixR()(i, i);
    min_r            = std::min(rii, min_r);
    max_r            = std::max(rii, max_r);
  }
  double rcond = min_r / max_r;
  if (rcond < 0 || std::isnan(rcond))
    rcond = 0;

  return rcond;
}

/// Deletes all dead directions from fullVector and returns a vector of reduced dimensionality.
inline double computeSquaredDifference(
    const std::array<double, 3> &u,
    std::array<double, 3>        v,
    const std::array<bool, 3>   &activeAxis = {{true, true, true}})
{
  // Subtract the values and multiply out dead dimensions
  for (unsigned int d = 0; d < v.size(); ++d) {
    v[d] = (u[d] - v[d]) * static_cast<int>(activeAxis[d]);
  }
  // @todo: this can be replaced by std::hypot when moving to C++17
  return std::accumulate(v.begin(), v.end(), static_cast<double>(0.), [](auto &res, auto &val) { return res + val * val; });
}

} // namespace precice::utils
