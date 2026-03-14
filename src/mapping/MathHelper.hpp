#pragma once

/**
 * @file MathHelper.hpp
 * @brief Low-level linear algebra helper: blockwise inversion of a lower triangular matrix.
 *
 * ## Why blockwise inversion?
 *
 * Given the Cholesky decomposition C = L * L^T of the RBF kernel matrix,
 * the diagonal of C^{-1} is needed for Rippa's Leave-One-Out Cross-Validation
 * (LOOCV) error estimate. To compute it we need the inverse of the lower
 * triangular factor L.
 *
 * A naive approach calls Eigen's `inverse()` which computes the full n×n
 * inverse matrix (O(n^3) time, O(n^2) memory). For large RBF systems this
 * is unnecessarily expensive; we only need the diagonal entries of C^{-1},
 * which are the column-squared-sum of L^{-1}.
 *
 * **Blockwise triangular solve** reduces the constant factor by exploiting
 * the triangular structure:
 *   - Diagonal blocks are inverted independently (small triangular solves).
 *   - Off-diagonal blocks are computed via a matrix-multiply recurrence.
 *   - Cache locality is improved when the block size matches the L1 cache line.
 *
 * Block size selection mirrors Eigen's internal heuristic:
 *   blockSize = round_down_to_multiple_of_16(n / 8), clamped to [8, 128].
 *
 * ## Block Inversion Recurrence
 *
 * Partition L into blocks. For the (i, j) block of L^{-1} (i > j):
 *
 *   L_inv(i, j) = -L_inv(i, i) * [ L(i,j)*L_inv(j,j) + Σ_{k=j+1}^{i-1} L(i,k)*L_inv(k,j) ]
 *
 * This follows from the standard block forward-substitution identity:
 *   L * L_inv = I  ⟹  L(i,i)*L_inv(i,j) + Σ_{k<i} L(i,k)*L_inv(k,j) = 0  (for i > j)
 *   ⟹  L_inv(i,j) = -L_inv(i,i) * Σ_{k=j}^{i-1} L(i,k)*L_inv(k,j)
 */

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

  // BLOCK SIZE SELECTION
  // Mirror Eigen's internal heuristic: blockSize = n/8, rounded down to the
  // nearest multiple of 16, then clamped to [8, 128]. This ensures the blocks
  // fit comfortably in L1/L2 cache for typical matrix dimensions encountered
  // in preCICE RBF problems.
  int blockSize = static_cast<int>(n / 8);
  blockSize     = (blockSize / 16) * 16;              // round down to multiple of 16
  blockSize     = std::min(std::max(blockSize, 8), 128); // clamp to [8, 128]

  // Edge case: if n < 8, blockSize would be 0 — treat the whole matrix as one block.
  if (blockSize == 0) {
    blockSize = n;
  }
  const int numBlocks = static_cast<int>((n + blockSize - 1) / blockSize); // ceil division

  // Initialize L_inv as the identity matrix.
  // This serves double duty: diagonal blocks start as "I" for in-place triangular solve,
  // and off-diagonal entries start at 0 (correct initial value for the recurrence).
  // TODO: maybe consider to return by reference, as Eigen doesn't recommend to return by value
  Eigen::MatrixXd L_inv = Eigen::MatrixXd::Identity(n, n);

  // Iterate over diagonal blocks from top-left to bottom-right.
  for (int i = 0; i < numBlocks; ++i) {
    int i_start     = i * blockSize;
    int i_end       = std::min(i_start + blockSize, n);
    int i_blockSize = i_end - i_start; // actual block size (may be smaller for the last block)

    // STEP 1: Invert the i-th diagonal block of L.
    // L(i_block, i_block) is a small lower-triangular matrix.
    // solveInPlace(I_block) computes: L_block^{-1} * I_block = L_block^{-1}
    // and stores the result back into L_inv(i_block, i_block).
    // Each diagonal block is inverted independently — O(blockSize^3) per block.
    L.block(i_start, i_start, i_blockSize, i_blockSize)
        .triangularView<Eigen::Lower>()
        .solveInPlace(L_inv.block(i_start, i_start, i_blockSize, i_blockSize));

    // STEP 2: Compute the off-diagonal blocks L_inv(i, j) for j < i.
    // Using the block forward-substitution recurrence:
    //   L_inv(i, j) = -L_inv(i,i) * [ L(i,j)*L_inv(j,j) + Σ_{k=j+1}^{i-1} L(i,k)*L_inv(k,j) ]
    //
    // We build this sum incrementally in `offDiagBlock`, starting with the
    // j-th term, then accumulating the intermediate k-terms in the inner loop.
    for (int j = 0; j < i; ++j) {
      int j_start     = j * blockSize;
      int j_end       = std::min(j_start + blockSize, n);
      int j_blockSize = j_end - j_start;

      // Term for k = j: L(i, j) * L_inv(j, j)
      // This is the "direct" sub-diagonal coupling from block j to block i.
      Eigen::MatrixXd offDiagBlock = L.block(i_start, j_start, i_blockSize, j_blockSize) * L_inv.block(j_start, j_start, j_blockSize, j_blockSize)
                                                                                               .triangularView<Eigen::Lower>();

      // Accumulate the sum of L(i, k)*L_inv(k, j)  for k = j+1 ... i-1.
      // These are "indirect" couplings passing through intermediate blocks.
      for (int k = j + 1; k < i; ++k) {
        int k_start     = k * blockSize;
        int k_end       = std::min(k_start + blockSize, n);
        int k_blockSize = k_end - k_start;

        // Add the L(i, k) * L_inv(k, j) contribution.
        // Note: L_inv(k, j) is already computed (we process blocks left-to-right).
        offDiagBlock += L.block(i_start, k_start, i_blockSize, k_blockSize) * L_inv.block(k_start, j_start, k_blockSize, j_blockSize);
      }

      // Final recurrence step:
      //   L_inv(i, j) = -L_inv(i, i) * (accumulated sum)
      // The negative sign comes from moving the off-diagonal L terms to the RHS of L*L_inv = 0.
      // Note: although L_inv appears on both left and right, the SOURCE (i_start rows of offDiagBlock)
      // and DESTINATION (i_block of L_inv at column j) are distinct — no aliasing issue.
      L_inv.block(i_start, j_start, i_blockSize, j_blockSize) = -L_inv.block(i_start, i_start, i_blockSize, i_blockSize) * offDiagBlock;
    }
  }
  return L_inv;
}
} // namespace precice::utils
