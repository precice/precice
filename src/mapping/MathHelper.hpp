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
inline Eigen::MatrixXd invertLowerTriangularBlockwise(const Eigen::MatrixXd &L)
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
  Eigen::MatrixXd L_inv = Eigen::MatrixXd::Identity(n, n);

  // Now we iterate over all blocks...
  for (int i = 0; i < numBlocks; ++i) {
    int i_start     = i * blockSize;
    int i_end       = std::min(i_start + blockSize, n);
    int i_blockSize = i_end - i_start;

    // Step 1: Invert the diagonal block
    L.block(i_start, i_start, i_blockSize, i_blockSize)
        .triangularView<Eigen::Lower>()
        .solveInPlace(L_inv.block(i_start, i_start, i_blockSize, i_blockSize));

    // Step 2: Compute off-diagonal blocks
    for (int j = 0; j < i; ++j) {
      int j_start     = j * blockSize;
      int j_end       = std::min(j_start + blockSize, n);
      int j_blockSize = j_end - j_start;

      // computation for the off-diagonal block using the block Lij and Ljj (from L_inv)
      Eigen::MatrixXd offDiagBlock = L.block(i_start, j_start, i_blockSize, j_blockSize) * L_inv.block(j_start, j_start, j_blockSize, j_blockSize)
                                                                                               .triangularView<Eigen::Lower>();
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
} // namespace precice::utils