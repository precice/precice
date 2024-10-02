#pragma once

#include <Eigen/Core>

namespace precice::utils {
inline Eigen::MatrixXd invertLowerTriangularBlockwise(const Eigen::MatrixXd& L) {
    int n = L.rows();
    int blockSize = 64;  // Block size can be adjusted based on cache size and performance

    // Initialize the inverse matrix as a zero matrix
    Eigen::MatrixXd L_inv = Eigen::MatrixXd::Zero(n, n);

    // Number of blocks
    int numBlocks = (n + blockSize - 1) / blockSize;

    // Store inverses of diagonal blocks
    std::vector<Eigen::MatrixXd> diagInverses(numBlocks);

    // Step 1: Invert diagonal blocks
    for (int k = 0; k < numBlocks; ++k) {
        int start = k * blockSize;
        int end = std::min(start + blockSize, n);
        int currentBlockSize = end - start;

        // Extract the diagonal block Lkk
        auto Lkk = L.block(start, start, currentBlockSize, currentBlockSize);

        // Invert the diagonal block
        Eigen::MatrixXd Lkk_inv = Lkk.triangularView<Eigen::Lower>().solve(
            Eigen::MatrixXd::Identity(currentBlockSize, currentBlockSize));

        // Place the inverted diagonal block into L_inv
        L_inv.block(start, start, currentBlockSize, currentBlockSize) = Lkk_inv;

        // Store the inverted diagonal block for later use
        diagInverses[k] = Lkk_inv;
    }
// Step 2: Compute off-diagonal blocks
for (int i = 1; i < numBlocks; ++i) {
    int i_start = i * blockSize;
    int i_end = std::min(i_start + blockSize, n);
    int i_blockSize = i_end - i_start;

    // Inverse of the diagonal block Lii (already computed)
    Eigen::MatrixXd Lii_inv = diagInverses[i];

    for (int j = 0; j < i; ++j) {
        int j_start = j * blockSize;
        int j_end = std::min(j_start + blockSize, n);
        int j_blockSize = j_end - j_start;

        // Initialize the off-diagonal block
        Eigen::MatrixXd offDiagBlock = Eigen::MatrixXd::Zero(i_blockSize, j_blockSize);

        // Accumulate the sum of L_ik * L_inv(k, j) for k = j+1 to i-1
        for (int k = j + 1; k < i; ++k) {
            int k_start = k * blockSize;
            int k_end = std::min(k_start + blockSize, n);
            int k_blockSize = k_end - k_start;

            // L_ik
            Eigen::MatrixXd Lik = L.block(i_start, k_start, i_blockSize, k_blockSize);

            // L_inv(k, j), previously computed
            Eigen::MatrixXd L_inv_kj = L_inv.block(k_start, j_start, k_blockSize, j_blockSize);

            // Update offDiagBlock
            offDiagBlock += Lik * L_inv_kj;
        }

        // Extract the block Lij
        Eigen::MatrixXd Lij = L.block(i_start, j_start, i_blockSize, j_blockSize);

        // Inverse of the diagonal block Ljj (already computed)
        Eigen::MatrixXd Ljj_inv = diagInverses[j];

        // Final computation for the off-diagonal block
        offDiagBlock += Lij * Ljj_inv;
        offDiagBlock = -Lii_inv * offDiagBlock;

        // Place the computed off-diagonal block into L_inv
        L_inv.block(i_start, j_start, i_blockSize, j_blockSize) = offDiagBlock;
    }
}

    return L_inv;
}

// inline void invertLowerTriangular(const Eigen::Ref<const Eigen::MatrixXd>& L, Eigen::MatrixXd& L_inv)
// {
//     int n = L.rows();
//     L_inv.resize(n, n);  // Ensure L_inv has the correct size

//     if (n <= 16) { // Base case threshold
//         // Use Eigen's dense inverse for small matrices
//         L_inv = L.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(n, n));
//     } else {
//         int k = n / 2;

//         // Initialize L_inv to zero
//         L_inv.setZero();

//         // Recursively invert diagonal blocks
//         invertLowerTriangular(L.topLeftCorner(k, k), L_inv.topLeftCorner(k, k));
//         invertLowerTriangular(L.bottomRightCorner(n - k, n - k), L_inv.bottomRightCorner(n - k, n - k));

//         // Compute the off-diagonal block
//         Eigen::MatrixXd L21 = L.bottomLeftCorner(n - k, k);
//         Eigen::MatrixXd temp = L_inv.bottomRightCorner(n - k, n - k) * L21;
//         L_inv.bottomLeftCorner(n - k, k) = -temp * L_inv.topLeftCorner(k, k);
//     }
// }

// template<typename MatrixType>
// MatrixType invertLowerTriangularUnblocked(const MatrixType& m)
// {
//     MatrixType invM = MatrixType::Identity(m.rows(), m.cols());
//     invM = m.template triangularView<Eigen::Lower>().solve(invM);
//     return invM;
// }

// template<typename MatrixType>
// MatrixType invertLowerTriangular(const MatrixType& m)
// {
//     typedef typename MatrixType::Index Index;
//     Index size = m.rows();


//     Index n1 = blockSize;
//     Index n2 = size - n1;

//     // Partition the matrix into blocks
//     Eigen::MatrixXd m11 = m.topLeftCorner(n1, n1);
//     Eigen::MatrixXd m21 = m.bottomLeftCorner(n2, n1);
//     Eigen::MatrixXd m22 = m.bottomRightCorner(n2, n2);

//     // Recursively compute the inverses of m11 and m22
//     auto invM11 = invertLowerTriangular(m11);
//     auto invM22 = invertLowerTriangular(m22);

//     // Compute invM21 = -invM22 * m21 * invM11
//     auto invM21 = -invM22 * m21 * invM11;

//     // Assemble the inverse matrix
//     MatrixType invM(size, size);
//     invM.setZero();

//     invM.topLeftCorner(n1, n1) = invM11;
//     invM.bottomLeftCorner(n2, n1) = invM21;
//     invM.bottomRightCorner(n2, n2) = invM22;

//     return invM;
// }

// inline Eigen::MatrixXd invertLowerTriangular(const Eigen::MatrixXd &L)
// {
//     int n = L.rows();

//     // Base case: if the matrix is small, invert directly
//     if (n < 32) {
//         // Use Eigen's dense inverse for small matrices
//         return L.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(n, n));
//     } else {
//         // Compute blockSize based on the matrix size and hardware considerations
//         int blockSize = n / 8;
//         blockSize = (blockSize / 16) * 16; // Align blockSize to multiple of 16
//         blockSize = std::min(std::max(blockSize, 8), 128); // Clamp blockSize between 8 and 128

//         // Ensure blockSize is valid
//         if (blockSize == 0 || blockSize >= n)
//             blockSize = n / 2;

//         int k = blockSize;
//         int n_k = n - k;

//         // Partition the matrix into blocks
//         Eigen::MatrixXd L11 = L.topLeftCorner(k, k);
//         Eigen::MatrixXd L21 = L.bottomLeftCorner(n_k, k);
//         Eigen::MatrixXd L22 = L.bottomRightCorner(n_k, n_k);

//         // Recursively invert diagonal blocks
//         Eigen::MatrixXd L11_inv = invertLowerTriangular(L11);
//         Eigen::MatrixXd L22_inv = invertLowerTriangular(L22);

//         // Compute the off-diagonal block
//         Eigen::MatrixXd L21_inv = -L22_inv * L21 * L11_inv;

//         // Assemble the inverse matrix
//         Eigen::MatrixXd L_inv(n, n);
//         L_inv.setZero();
//         L_inv.topLeftCorner(k, k)             = L11_inv;
//         L_inv.bottomLeftCorner(n_k, k)        = L21_inv;
//         L_inv.bottomRightCorner(n_k, n_k)     = L22_inv;

//         return L_inv;
//     }
// }
} // namespace precice::utils