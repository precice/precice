#pragma once

/**
 * @file RadialBasisFctSolver.hpp
 * @brief Dense CPU solver for Radial Basis Function (RBF) interpolation.
 *
 * ## Purpose
 * `RadialBasisFctSolver<RBF_T>` assembles and caches the RBF interpolation matrices
 * during mapping setup (once, O(n^3)) and performs fast backward-substitution
 * during each data exchange (reuse, O(n^2)).
 *
 * ## Workflow
 * 1. **Setup (Constructor)**:
 *    - Assemble the symmetric kernel matrix C: C_ij = phi(||x_i - x_j||)
 *    - Decompose C: Cholesky for SPD kernels, QR for general kernels
 *    - Assemble the evaluation matrix A: A_ij = phi(||y_i - x_j||)
 *    - Optionally assemble polynomial matrices Q (input) and V (output)
 *
 * 2. **Mapping (solveConsistent / solveConservative)**:
 *    - Consistent:   compute lambda = C^{-1} * data, then output = A * lambda
 *    - Conservative: compute Au = A^T * data, then output = C^{-T} * Au
 *
 * ## Just-In-Time (JIT) Caching
 * For on-the-fly mapping (when the receiving mesh isn't known at setup time),
 * the solver exposes `computeCacheData`, `interpolateAt`, `addWriteDataToCache`,
 * and `evaluateConservativeCache`. These split the solve into a precomputation step
 * (independent of query locations) and a per-query evaluation step, allowing
 * efficient repeated evaluation at arbitrary points after a single matrix solve.
 */

#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/irange.hpp>
#include <numeric>
#include <type_traits>
#include "mapping/MathHelper.hpp"
#include "mapping/config/MappingConfigurationTypes.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"

namespace precice::mapping {

/**
 * This class assembles and solves an RBF system, given an input mesh and an output mesh with relevant vertex IDs.
 * The class uses a dense matrix decomposition in order to decompose the resulting system(s) and a backward substitution
 * in order to solve the system at runtime. The functionality uses Eigen and supports only serial execution. In case
 * the polynomial="separate" option is used, the polynomial system is solved using a QR decomposition.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class RadialBasisFctSolver {
public:
  using DecompositionType = std::conditional_t<RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite(), Eigen::LLT<Eigen::MatrixXd>, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>;
  using BASIS_FUNCTION_T  = RADIAL_BASIS_FUNCTION_T;
  /// Default constructor
  RadialBasisFctSolver() = default;

  /**
   * assembles the system matrices and computes the decomposition of the interpolation matrix
   * inputMesh refers to the mesh where the interpolants are built on, i.e., the input mesh
   * for consistent mappings and the output mesh for conservative mappings
   * outputMesh refers to the mesh where we evaluate the interpolants, i.e., the output mesh
   * consistent mappings and the input mesh for conservative mappings
   */
  template <typename IndexContainer>
  RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                       const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial);

  /// Maps the given input data
  Eigen::MatrixXd solveConsistent(Eigen::MatrixXd &inputData, Polynomial polynomial) const;

  void computeCacheData(Eigen::MatrixXd &inputData, Polynomial polynomial, Eigen::MatrixXd &polyOut, Eigen::MatrixXd &coeffsOut) const;

  /// Maps the given input data
  Eigen::MatrixXd solveConservative(const Eigen::MatrixXd &inputData, Polynomial polynomial) const;

  // Clear all stored matrices
  void clear();

  // Returns the size of the input data
  Eigen::Index getInputSize() const;

  // Returns the size of the input data
  Eigen::Index getOutputSize() const;

  // Returns the number of active axis in use
  Eigen::Index getNumberOfPolynomials() const;

  template <typename IndexContainer>
  Eigen::VectorXd interpolateAt(const mesh::Vertex &v, const Eigen::MatrixXd &poly, const Eigen::MatrixXd &coeffs,
                                const RADIAL_BASIS_FUNCTION_T &function, const IndexContainer &inputIDs, const mesh::Mesh &inMesh) const;

  template <typename IndexContainer>
  void addWriteDataToCache(const mesh::Vertex &v, const Eigen::VectorXd &load, Eigen::MatrixXd &epsilon, Eigen::MatrixXd &Au,
                           const RADIAL_BASIS_FUNCTION_T &basisFunction, const IndexContainer &inputIDs, const mesh::Mesh &inMesh) const;

  void evaluateConservativeCache(Eigen::MatrixXd &epsilon, const Eigen::MatrixXd &Au, Eigen::MatrixXd &result) const;

private:
  mutable precice::logging::Logger _log{"mapping::RadialBasisFctSolver"};

  double evaluateRippaLOOCVerror(const Eigen::VectorXd &lambda) const;
  /// Decomposition of the interpolation matrix
  DecompositionType _decMatrixC;

  /// Diagonal entris of the inverse matrix C, requires for the Rippa scheme
  Eigen::VectorXd _inverseDiagonal;

  /// Decomposition of the polynomial (for separate polynomial)
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _qrMatrixQ;

  /// Polynomial matrix of the input mesh (for separate polynomial)
  Eigen::MatrixXd _matrixQ;

  /// Polynomial matrix of the output mesh (for separate polynomial)
  Eigen::MatrixXd _matrixV;

  /// Evaluation matrix (output x input)
  Eigen::MatrixXd _matrixA;

  bool                computeCrossValidation = false;
  std::array<bool, 3> _localActiveAxis;
};

// ------- Non-Member Functions ---------

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
  // use std::inner_product to avoid the need for lambda expressions
  return std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
}

/// given the active axis, computes sets the axis with the lowest spatial expansion to dead
template <typename IndexContainer>
constexpr void reduceActiveAxis(const mesh::Mesh &mesh, const IndexContainer &IDs, std::array<bool, 3> &axis)
{
  // make a pair of the axis and the difference
  std::array<std::pair<int, double>, 3> differences;

  // Compute the difference magnitude per direction
  for (std::size_t d = 0; d < axis.size(); ++d) {
    // Ignore dead axis here, i.e., apply the max value such that they are sorted on the last position(s)
    if (axis[d] == false) {
      differences[d] = std::make_pair<int, double>(d, std::numeric_limits<double>::max());
    } else {
      auto res = std::minmax_element(IDs.begin(), IDs.end(), [&](const auto &a, const auto &b) { return mesh.vertex(a).coord(d) < mesh.vertex(b).coord(d); });
      // Check if we are above or below the threshold
      differences[d] = std::make_pair<int, double>(d, std::abs(mesh.vertex(*res.second).coord(d) - mesh.vertex(*res.first).coord(d)));
    }
  }

  std::sort(differences.begin(), differences.end(), [](const auto &d1, const auto &d2) { return d1.second < d2.second; });
  // Disable the axis having the smallest expansion
  axis[differences[0].first] = false;
}

// Fill in the polynomial entries
template <typename IndexContainer>
inline void fillPolynomialEntries(Eigen::MatrixXd &matrix, const mesh::Mesh &mesh, const IndexContainer &IDs, Eigen::Index startIndex, std::array<bool, 3> activeAxis)
{
  // Loop over all vertices in the mesh
  for (const auto &i : IDs | boost::adaptors::indexed()) {

    // 1. the constant contribution
    matrix(i.index(), startIndex) = 1.0;

    // 2. the linear contribution
    const auto  &u = mesh.vertex(i.value()).rawCoords();
    unsigned int k = 0;
    // Loop over all three space dimension and ignore dead axis
    for (unsigned int d = 0; d < activeAxis.size(); ++d) {
      if (activeAxis[d]) {
        PRECICE_ASSERT(matrix.rows() > i.index(), matrix.rows(), i.index());
        PRECICE_ASSERT(matrix.cols() > startIndex + 1 + k, matrix.cols(), startIndex + 1 + k);
        matrix(i.index(), startIndex + 1 + k) = u[d];
        ++k;
      }
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                               std::array<bool, 3> activeAxis, Polynomial polynomial)
{
  /**
   * BUILD RBF INTERPOLATION MATRIX (Matrix C or CLU)
   *
   * Constructs the core matrix for RBF interpolation. This is the Gram matrix of the
   * radial basis function evaluated at all input data points.
   *
   * Mathematical formulation:
   * For strictly positive definite functions:
   *   C_ij = phi(||x_i - x_j||)  for all input vertices i,j
   *   This matrix is symmetric (C = C^T) and positive definite.
   *
   * If polynomial="on":
   *   Extended system augments C with polynomial basis:
   *   | C  P |
   *   | P^T 0 |
   *   where P contains polynomial evaluations (1, x, y, z for each vertex)
   *
   * If polynomial="separate":
   *   Only builds C; polynomial handled separately in dedicated matrices Q and V.
   *
   * Dead axis handling:
   * For 2D problems, one axis is "dead" (all vertices have same coordinate).
   * Distance computation ignores dead axes to avoid numerical issues.
   */

  // Treat the 2D case as 3D case with dead axis
  // This simplifies logic: compute in 3D but ignore zero-extent dimensions
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  const unsigned int dimensions     = 3;
  const unsigned int polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  // Add linear polynom degrees if polynomial requires this
  const auto inputSize = inputIDs.size();
  const auto n         = inputSize + polyparams;

  PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  Eigen::MatrixXd matrixCLU(n, n);

  // Required to fill the poly -> poly entries in the matrix, which remain otherwise untouched
  if (polynomial == Polynomial::ON) {
    matrixCLU.setZero();
  }

  /**
   * COMPUTE RBF KERNEL MATRIX ENTRIES
   *
   * Evaluate the radial basis function phi for all pairs of input vertices:
   * C_ij = phi(r_ij) where r_ij = ||x_i - x_j|| (Euclidean distance)
   *
   * Kernel examples:
   * - Multiquadric: phi(r) = sqrt(r^2 + c^2)
   * - Inverse multiquadric: phi(r) = 1 / sqrt(r^2 + c^2)
   * - Thin-plate spline: phi(r) = r^2 log(r)
   * - Gaussian: phi(r) = exp(-c*r^2)
   *
   * Symmetry optimization: Only compute upper triangle, copy to lower.
   */
  auto         i_iter  = inputIDs.begin();
  Eigen::Index i_index = 0;
  for (; i_iter != inputIDs.end(); ++i_iter, ++i_index) {
    const auto &u       = inputMesh.vertex(*i_iter).rawCoords();
    auto        j_iter  = i_iter;
    auto        j_index = i_index;
    for (; j_iter != inputIDs.end(); ++j_iter, ++j_index) {
      const auto &v                 = inputMesh.vertex(*j_iter).rawCoords();
      // Compute squared distance, respecting dead axes (ignore dimensions with zero extent)
      double      squaredDifference = computeSquaredDifference(u, v, activeAxis);
      // Evaluate kernel at this distance
      matrixCLU(i_index, j_index)   = basisFunction.evaluate(std::sqrt(squaredDifference));
    }
  }

  /**
   * ADD POLYNOMIAL BASIS ROWS/COLUMNS (if integrated polynomial is enabled)
   *
   * The polynomial part adds linear degrees of freedom to the RBF system.
   * This allows the interpolant to exactly reproduce linear functions.
   *
   * Polynomial contributions:
   * P = [p_0, p_1, ..., p_k]  for each vertex, where:
   *   p_0(x) = 1 (constant)
   *   p_1(x) = x_1 (linear in first active dimension)
   *   p_2(x) = x_2 (linear in second active dimension)
   *   ... (only for active axes, dead axes are skipped)
   *
   * Block structure of extended matrix C_ext:
   * | RBF_part   Poly_part |   <- RBF evaluations, then polynomial columns
   * | Poly_part^T    0     |   <- Polynomial rows, zero in lower-right (singular submatrix)
   */
  // Add potentially the polynomial contribution in the matrix
  if (polynomial == Polynomial::ON) {
    fillPolynomialEntries(matrixCLU, inputMesh, inputIDs, inputSize, activeAxis);
  }
  // Symmetry: copy upper triangle to lower triangle (matrix is symmetric)
  matrixCLU.triangularView<Eigen::Lower>() = matrixCLU.transpose();
  return matrixCLU;
}

template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
Eigen::MatrixXd buildMatrixA(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                             const mesh::Mesh &outputMesh, const IndexContainer outputIDs, std::array<bool, 3> activeAxis, Polynomial polynomial)
{
  /**
   * BUILD RBF EVALUATION MATRIX (Matrix A)
   *
   * This matrix evaluates the RBF basis functions at output locations using input data.
   * Dimensions: (num_output_vertices) x (num_input_vertices + num_polynomial_params)
   *
   * Mathematical role:
   * Once we solve for RBF coefficients lambda from: C * [lambda; poly_coeff] = [data; 0]
   * We evaluate the RBF at output points using: A * [lambda; poly_coeff] = output
   *
   * A_ij = phi(||y_i - x_j||)  for output vertex i and input vertex j
   * where phi is the radial basis function (same as used in matrix C).
   */

  // Treat the 2D case as 3D case with dead axis
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  const unsigned int dimensions     = 3;
  const unsigned int polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  const auto inputSize  = inputIDs.size();
  const auto outputSize = outputIDs.size();
  const auto n          = inputSize + polyparams;

  PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  Eigen::MatrixXd matrixA(outputSize, n);

  /**
   * COMPUTE RBF VALUES FOR ALL OUTPUT-INPUT PAIRS
   *
   * For each output vertex, evaluate the basis function at each input vertex.
   * A_ij = phi(distance between output_i and input_j)
   */
  for (const auto &i : outputIDs | boost::adaptors::indexed()) {
    const auto &u = outputMesh.vertex(i.value()).rawCoords();
    for (const auto &j : inputIDs | boost::adaptors::indexed()) {
      const auto &v                 = inputMesh.vertex(j.value()).rawCoords();
      // Compute squared distance, respecting dead axes
      double      squaredDifference = computeSquaredDifference(u, v, activeAxis);
      matrixA(i.index(), j.index()) = basisFunction.evaluate(std::sqrt(squaredDifference));
    }
  }

  /**
   * ADD POLYNOMIAL BASIS COLUMNS (if integrated polynomial is enabled)
   *
   * Similar to matrix C, we add polynomial columns to enable exact reproduction
   * of linear functions at output locations.
   */
  // Add potentially the polynomial contribution in the matrix
  if (polynomial == Polynomial::ON) {
    fillPolynomialEntries(matrixA, outputMesh, outputIDs, inputSize, activeAxis);
  }
  return matrixA;
}

// Variant operating on the Cholesky decomposition
inline Eigen::VectorXd computeInverseDiagonal(Eigen::LLT<Eigen::MatrixXd> decMatrixC)
{
  /**
   * EFFICIENT COMPUTATION OF INVERSE MATRIX DIAGONAL (Cholesky variant)
   *
   * For Leave-One-Out Cross Validation (LOOCV) error estimation, we need the
   * diagonal of C^{-1}. Computing the full inverse is O(n^3); this method
   * computes only the diagonal in O(n^2) using the Cholesky decomposition.
   *
   * Mathematical approach:
   * If C = L*L^T (Cholesky decomposition), then:
   *   C^{-1} = (L^T)^{-1} * L^{-1}
   *
   * Diagonal entries:
   *   [C^{-1}]_ii = sum_k [L^{-T}]_ik * [L^{-1}]_ki = sum_k [L^{-1}]_ki^2
   *
   * Algorithm:
   * 1. Compute L^{-1} by solving L*X = I (blockwise for efficiency)
   * 2. Sum squares of each column to get diagonal entries
   *
   * Complexity: O(n^2) vs O(n^3) for full inverse
   * Benefit: Only 4 values per LOOCV computation needed.
   */

  // Compute the inverse of the lower triangular matrix L
  // We use a blockwise triangular inverse (more efficient than dense inverse)
  // This solves L * L_inv = I using an optimized triangular solution algorithm
  Eigen::MatrixXd L_inv = utils::invertLowerTriangularBlockwise(decMatrixC.matrixL());

  // Compute the diagonal elements of A^{-1} by evaluating (L^T)^{-1}L^{-1}
  // Each diagonal entry is the sum of squares along the corresponding column
  Eigen::VectorXd inverseDiagonal = (L_inv.array().square().colwise().sum()).transpose();

  return inverseDiagonal;
}

// Variant operating on the QR decomposition
inline Eigen::VectorXd computeInverseDiagonal(Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decMatrixC)
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

template <typename RADIAL_BASIS_FUNCTION_T>
double RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::evaluateRippaLOOCVerror(const Eigen::VectorXd &lambda) const
{
  // Implementation of LOOCV according to Rippa(1999), DOI: 10.1023/a:1018975909870
  double             loocv = 0;
  const Eigen::Index n     = lambda.size();
  // 2: Next, compute the RBF coefficient. These are exactly the same as computed during
  // the solve_consistent and we could also pass them into the function, depending on how
  // often wen want to compute the LOOCV. We omit the polynomial contribution here, i.e.,
  // the input data remains untouched
  // Eigen::VectorXd lambda = _decMatrixC.solve(inputData);

  // 3: Evaluate the Rippa formula:
  // The error estimate is given by a component-wise division: lambda/A^{-1}_{ii}.
  // We then compute the RMS of all LOOCV error entries (other options for the
  // aggregation should be possible)
  // TODO:Consider storing the reciprocal inverse diagonal entries
  loocv = std::sqrt((lambda.array() / _inverseDiagonal.array()).array().square().sum() / n);

  return loocv;
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                                                                    const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial)
{
  /**
   * RBF SOLVER CONSTRUCTOR - SYSTEM ASSEMBLY AND DECOMPOSITION
   *
   * This constructor is the core of RBF setup. It:
   * 1. Assembles the interpolation matrix C from input data points
   * 2. Decomposes C into a form efficient for solving
   * 3. Assembles the evaluation matrix A for later interpolation
   * 4. Handles polynomial terms (if enabled)
   * 5. Validates numerical properties (invertibility, condition numbers)
   *
   * Key design:
   * - Matrix decompositions are performed once in constructor
   * - During mapping, we only perform fast backward-substitution (O(n^2))
   * - This amortizes expensive O(n^3) decomposition cost over many evaluations
   */

  PRECICE_ASSERT(!(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite() && polynomial == Polynomial::ON),
                 "The integrated polynomial (polynomial=\"on\") is not supported for the selected radial-basis function. "
                 "Please select another radial-basis function or change the polynomial configuration.");

  // Convert dead axis vector into an active axis array for cleaner indexing
  // activeAxis[i] = true means dimension i has varying coordinates (active)
  // activeAxis[i] = false means dimension i is constant across all vertices (dead)
  std::array<bool, 3> activeAxis({{false, false, false}});
  std::transform(deadAxis.begin(), deadAxis.end(), activeAxis.begin(),
    [](const auto ax) { return !ax; });

  /**
   * STEP 1: ASSEMBLE AND DECOMPOSE INTERPOLATION MATRIX
   *
   * Build the RBF kernel matrix C from input vertices and choose appropriate decomposition:
   *
   * For strictly positive definite basis functions:
   *   - Use Cholesky decomposition (C = L*L^T)
   *   - Assumes C is symmetric positive definite
   *   - Fast: O(n^3) but with small constants
   *   - Solves: L*L^T*lambda = data (two triangular solves)
   *
   * For other basis functions (e.g., thin-plate spline):
   *   - Use QR decomposition with column pivoting
   *   - Handles rank-deficient or singular matrices
   *   - More robust to ill-conditioning
   *   - Used when C is not guaranteed positive definite
   */
  bool decompositionSuccessful = false;
  if constexpr (RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite()) {
    // Cholesky decomposition for positive definite kernels
    _decMatrixC             = buildMatrixCLU(basisFunction, inputMesh, inputIDs, activeAxis, polynomial).llt();
    decompositionSuccessful = _decMatrixC.info() == Eigen::ComputationInfo::Success;
  } else {
    // QR decomposition with column pivoting for general kernels
    _decMatrixC             = buildMatrixCLU(basisFunction, inputMesh, inputIDs, activeAxis, polynomial).colPivHouseholderQr();
    decompositionSuccessful = _decMatrixC.isInvertible();
  }

  PRECICE_CHECK(decompositionSuccessful,
                "The interpolation matrix of the RBF mapping from mesh \"{}\" to mesh \"{}\" is not invertable. "
                "This means that the mapping problem is not well-posed. "
                "Please check if your coupling meshes are correct (e.g. no vertices are duplicated) or reconfigure "
                "your basis-function (e.g. reduce the support-radius).",
                inputMesh.getName(), outputMesh.getName());

  /**
   * OPTIONAL STEP: COMPUTE INVERSE DIAGONAL FOR CROSS-VALIDATION
   *
   * Leave-One-Out Cross Validation (LOOCV) estimates RBF quality without test data.
   * Rippa(1999) formula allows efficient LOOCV using only the inverse diagonal:
   *   error_i = lambda_i / [C^{-1}]_ii
   *
   * Disabled for integrated polynomials (system size becomes ambiguous).
   */
  if (polynomial != Polynomial::ON && computeCrossValidation) {
    precice::profiling::Event e("map.rbf.computeLOOCV");
    _inverseDiagonal = computeInverseDiagonal(_decMatrixC);
  }

  /**
   * STEP 2: ASSEMBLE EVALUATION MATRIX
   *
   * Build matrix A that maps RBF coefficients to output vertex values.
   * This is computed once; during mapping we evaluate A*lambda efficiently.
   */
  _matrixA = buildMatrixA(basisFunction, inputMesh, inputIDs, outputMesh, outputIDs, activeAxis, polynomial);

  /**
   * STEP 3: HANDLE POLYNOMIAL TERMS (if polynomial="separate")
   *
   * For separated polynomial handling:
   * - Build matrix Q from polynomial basis evaluated at input vertices
   * - Build matrix V from polynomial basis evaluated at output vertices
   * - Decompose Q using QR for efficient polynomial solving
   *
   * Separated vs integrated polynomial:
   * - Integrated: polynomial added to RBF system (augmented matrix)
   * - Separated: polynomial solved independently, subtracted from data before RBF
   *            then added back to result after RBF evaluation
   *
   * Advantage of separated: Avoids ill-conditioning of augmented matrix blocks
   * (RBF block and polynomial block have very different scales)
   *
   * Condition number check:
   * If polynomial basis has high condition number (> 1e5), one axis is disabled
   * (e.g., in 2D, use constant + x term, drop y term for degenerate cases)
   */
  // In case we deal with separated polynomials, we need dedicated matrices for the polynomial contribution
  if (polynomial == Polynomial::SEPARATE) {
    // 4 = 1 + dimensions(3) = maximum number of polynomial parameters
    _localActiveAxis        = activeAxis;
    unsigned int polyParams = getNumberOfPolynomials();

    // Iteratively check polynomial conditioning, disabling axes as needed
    do {
      // Build polynomial matrix Q for input vertices
      _matrixQ.resize(inputIDs.size(), polyParams);
      fillPolynomialEntries(_matrixQ, inputMesh, inputIDs, 0, _localActiveAxis);

      // Check condition number: ratio of max to min singular values
      // High condition number indicates numerical instability
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(_matrixQ);
      PRECICE_ASSERT(svd.singularValues().size() > 0);
      PRECICE_DEBUG("Singular values in polynomial solver: {}", svd.singularValues());
      const double conditionNumber = svd.singularValues()(0) /
                                     std::max(svd.singularValues()(svd.singularValues().size() - 1), math::NUMERICAL_ZERO_DIFFERENCE);
      PRECICE_DEBUG("Condition number: {}", conditionNumber);

      // If condition number too high, disable the axis with smallest span
      if (conditionNumber > 1e5) {
        reduceActiveAxis(inputMesh, inputIDs, _localActiveAxis);
        polyParams = getNumberOfPolynomials();
      } else {
        break;  // Condition number acceptable
      }
    } while (true);

    // Build polynomial matrix V for output vertices
    _matrixV.resize(outputIDs.size(), polyParams);
    fillPolynomialEntries(_matrixV, outputMesh, outputIDs, 0, _localActiveAxis);

    // Decompose Q for later polynomial solving
    _qrMatrixQ = _matrixQ.colPivHouseholderQr();
  }
}

/**
 * CONSERVATIVE SOLVE
 *
 * Maps data *conservatively* from output vertices (input data) to input vertices (output data).
 * Unlike consistent, conservative preserves the global integral (or sum) of the field.
 *
 * Mathematical derivation:
 *   Consistent mapping:    output = A * C^{-1} * input
 *   Conservative mapping:  output = (A * C^{-1})^T * input  = C^{-T} * A^T * input
 *   (since C is symmetric: C^{-T} = C^{-1})
 *
 * Algorithm:
 *   Step 1 (Polynomial, SEPARATE only):
 *     epsilon = V^T * inputData
 *     tau     = epsilon - Q^T * mu   (after step 2)
 *     sigma   = C^{-T} * tau        (correction term)
 *
 *   Step 2 (Core conservative solve):
 *     Au  = A^T * inputData        (project input onto RBF basis)
 *     mu  = C^{-1} * Au           (RBF coefficients for conservative mapping)
 *
 *   Step 3 (Correction, SEPARATE only):
 *     out = mu + sigma             (add polynomial correction)
 *
 * @note inputData must be a column-major matrix with rows = getOutputSize() (N_out vertices)
 * @param polynomial Polynomial handling mode (ON, OFF, SEPARATE)
 * @return Eigen::MatrixXd  Mapped values at input vertices
 */
template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::MatrixXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConservative(const Eigen::MatrixXd &inputData, Polynomial polynomial) const
{
  PRECICE_ASSERT((_matrixV.size() > 0 && polynomial == Polynomial::SEPARATE) || _matrixV.size() == 0, _matrixV.size());
  PRECICE_ASSERT(inputData.rows() == _matrixA.rows());

  // STEP 2a: Au = A^T * inputData
  //   This projects the output-mesh values back onto the RBF basis coefficients.
  //   Au has shape (inputSize + polyparams) x dataDims.
  Eigen::MatrixXd Au = _matrixA.transpose() * inputData;
  PRECICE_ASSERT(Au.rows() == _matrixA.cols());

  // STEP 2b: mu = C^{-1} * Au  (RBF coefficients)
  Eigen::MatrixXd out = _decMatrixC.solve(Au);

  if (polynomial == Polynomial::SEPARATE) {
    // STEP 3 (SEPARATE polynomial correction):
    // epsilon = V^T * inputData  — polynomial moments of the input values
    Eigen::MatrixXd epsilon = _matrixV.transpose() * inputData;
    PRECICE_ASSERT(epsilon.rows() == _matrixV.cols());

    // tau = epsilon - Q^T * mu  — polynomial residual after RBF projection
    epsilon -= _matrixQ.transpose() * out;
    PRECICE_ASSERT(epsilon.rows() == _matrixQ.cols());

    // sigma = C_poly^{-T} * tau  — correction via transposed polynomial QR
    // Subtracted (note the negative sign) to yield the final conservative coefficients
    out -= static_cast<Eigen::MatrixXd>(_qrMatrixQ.transpose().solve(-epsilon));
  }
  return out;
}

/**
 * CONSISTENT SOLVE
 *
 * Interpolates data from input to output vertices using the precomputed decomposition.
 *
 * Two-phase algorithm:
 *
 * Phase 1 (Polynomial subtraction — only if polynomial=SEPARATE):
 *   - Solve Q * poly_coeffs = inputData  using the precomputed QR of Q
 *   - Subtract the polynomial fit from input:  inputData -= Q * poly_coeffs
 *   - This removes the polynomial trend from the data before RBF interpolation.
 *
 * Phase 2 (RBF solve + evaluate):
 *   - Solve C * lambda = inputData  using the precomputed Cholesky or QR of C
 *   - Evaluate at output locations: output = A * lambda
 *
 * Phase 3 (Polynomial restoration — only if polynomial=SEPARATE):
 *   - Re-add the polynomial contribution at output locations: output += V * poly_coeffs
 *
 * @note inputData is modified in-place (polynomial subtracted if SEPARATE).
 *       The caller must not rely on its value after this call.
 * @param polynomial Polynomial handling mode (ON, OFF, SEPARATE)
 * @return Eigen::MatrixXd  Mapped values at output vertices; rows = outputSize
 */
template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::MatrixXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(Eigen::MatrixXd &inputData, Polynomial polynomial) const
{
  PRECICE_ASSERT((_matrixQ.size() > 0 && polynomial == Polynomial::SEPARATE) || _matrixQ.size() == 0);
  Eigen::MatrixXd polynomialContribution;

  // PHASE 1 (polynomial=SEPARATE only):
  // Factor out the polynomial part so the RBF system only sees the residual.
  // poly_coeffs = Q^{+} * inputData  (least-squares: best polynomial fit to input data)
  // residual    = inputData - Q * poly_coeffs
  if (polynomial == Polynomial::SEPARATE) {
    polynomialContribution = _qrMatrixQ.solve(inputData);
    inputData -= (_matrixQ * polynomialContribution);
  }

  // PHASE 2: Solve C * lambda = inputData (residual after polynomial removal)
  //   integrated polynomial (polynomial=ON): data already includes polynomial columns
  //   no polynomial (polynomial=OFF): straightforward RBF-only solve
  PRECICE_ASSERT(inputData.rows() == _matrixA.cols());
  Eigen::MatrixXd p = _decMatrixC.solve(inputData);

  // LOOCV cross-validation check (if enabled). Logs a warning if the mapping
  // is of poor quality based on the Rippa leave-one-out error estimate.
  if (polynomial != Polynomial::ON && computeCrossValidation) {
    precice::profiling::Event e("map.rbf.evaluateLOOCV");
    PRECICE_INFO("Cross validation error (LOOCV): {}", evaluateRippaLOOCVerror(p));
  }

  // Evaluate at output locations: output = A * lambda
  PRECICE_ASSERT(p.rows() == _matrixA.cols());
  Eigen::MatrixXd out = _matrixA * p;

  // PHASE 3 (polynomial=SEPARATE only):
  // Restore the polynomial contribution at output locations using matrix V:
  //   output += V * poly_coeffs
  if (polynomial == Polynomial::SEPARATE) {
    out += (_matrixV * polynomialContribution);
  }
  return out;
}

/**
 * COMPUTE CACHE DATA (JIT precomputation)
 *
 * Splits the consistent solve into a reusable precomputation step and a
 * per-query evaluation step. Called once per time step; the results (polyOut, coeffsOut)
 * are reused for multiple `interpolateAt()` calls at different query locations.
 *
 * This is the JIT equivalent of `solveConsistent`, but without the final A*lambda step.
 *
 * @param[in,out] inputData  Input values at input vertices (modified in-place; polynomial subtracted)
 * @param[in]     polynomial Polynomial handling mode
 * @param[out]    polyOut    Polynomial coefficients (meaningful only if SEPARATE)
 * @param[out]    coeffsOut  RBF lambda coefficients = C^{-1} * (inputData - Q * polyOut)
 */
template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::computeCacheData(Eigen::MatrixXd &inputData, Polynomial polynomial, Eigen::MatrixXd &polyOut, Eigen::MatrixXd &coeffsOut) const
{
  PRECICE_ASSERT((_matrixQ.size() > 0 && polynomial == Polynomial::SEPARATE) || _matrixQ.size() == 0);

  // Step 1 (polynomial=SEPARATE only): factor out polynomial part from the input.
  // polyOut = Q^{+} * inputData, then subtract Q * polyOut from input (residual).
  if (polynomial == Polynomial::SEPARATE) {
    polyOut = _qrMatrixQ.solve(inputData);
    inputData -= (_matrixQ * polyOut);
  }

  // Step 2: Solve C * coeffsOut = inputData  (precomputed decomposition, O(n^2))
  // coeffsOut holds the RBF lambda coefficients for the current input field.
  coeffsOut = _decMatrixC.solve(inputData);
}

/**
 * INTERPOLATE AT (JIT per-query evaluation — consistent)
 *
 * Given precomputed RBF coefficients (`coeffs`) and optional polynomial
 * coefficients (`poly`), evaluate the RBF interpolant at a single query vertex `v`.
 *
 * Formula:
 *   result = sum_j coeffs[j] * phi(||v - x_j||)     [RBF part]
 *          + poly[0]                                 [constant polynomial term]
 *          + sum_d poly[1+k] * v.coord[d]            [linear polynomial terms]
 *
 * This is the per-query evaluation partner to `computeCacheData`. Call
 * `computeCacheData` once per time step, then call `interpolateAt` for each query.
 *
 * @note Dead axes are NOT applied here because vertex coordinates are already
 *  in the correct 3D representation and no user-facing dead-axis option exists
 *  inside the PartitionOfUnity solver.
 *
 * @param v          Query vertex (location where we want the interpolated value)
 * @param poly       Polynomial coefficients from computeCacheData (empty if no SEPARATE poly)
 * @param coeffs     RBF lambda coefficients from computeCacheData
 * @param function   The RBF kernel (same as used in the constructor)
 * @param inputIDs   IDs of input vertices in `inMesh` (defines column ordering of coeffs)
 * @param inMesh     The input mesh containing the RBF source vertices
 * @return Eigen::VectorXd  Interpolated values at query vertex v (length = dataDims)
 */
template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
Eigen::VectorXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::interpolateAt(const mesh::Vertex &v, const Eigen::MatrixXd &poly, const Eigen::MatrixXd &coeffs,
                                                                              const RADIAL_BASIS_FUNCTION_T &basisFunction, const IndexContainer &inputIDs, const mesh::Mesh &inMesh) const
{
  PRECICE_TRACE();

  Eigen::VectorXd result(coeffs.cols());
  result.setZero();

  // RBF part: evaluate each input vertex's RBF contribution at query location v.
  // result += phi(||v - x_j||) * lambda_j  for each input vertex j.
  const auto &out = v.rawCoords();
  for (const auto &j : inputIDs | boost::adaptors::indexed()) {
    const auto &in                = inMesh.vertex(j.value()).rawCoords();
    double      squaredDifference = computeSquaredDifference(out, in);
    double      eval              = basisFunction.evaluate(std::sqrt(squaredDifference));
    result += eval * coeffs.row(j.index()).transpose();
  }

  // Polynomial part (SEPARATE mode only): add polynomial contribution at query point.
  //   poly.row(0): constant term coefficient  → multiplied by 1
  //   poly.row(k): linear term coefficient in dimension d → multiplied by v.coord[d]
  if (poly.size() > 0) {
    // constant polynomial
    result += 1 * poly.row(0).transpose();
    // the linear contributions — one per active dimension
    int k = 1;
    for (std::size_t d = 0; d < _localActiveAxis.size(); ++d) {
      if (_localActiveAxis[d]) {
        result += out[d] * poly.row(k).transpose();
        ++k;
      }
    }
  }
  return result;
}

/**
 * ADD WRITE DATA TO CACHE (JIT precomputation — conservative)
 *
 * Accumulates the conservative contribution of a single query vertex `v`
 * with load vector `load` into the running cache buffers `Au` and `epsilon`.
 *
 * This is the per-query step of the JIT conservative mapping. Call this for
 * every write vertex, then call `evaluateConservativeCache` once at the end.
 *
 * Mathematical role:
 *   Au_j      += phi(||v - x_j||) * load   for each input vertex j   [RBF column of A^T]
 *   epsilon_0 += 1 * load                                              [constant poly term]
 *   epsilon_k += v.coord[d] * load                                    [linear poly terms]
 *
 * @param v            Query vertex at which write data was prescribed
 * @param load         Data values at query vertex (scalar or vector per component)
 * @param[out] epsilon Accumulated polynomial moment matrix (rows = polyParams, cols = dataDims)
 * @param[out] Au      Accumulated A^T * load matrix (rows = inputIDs.size(), cols = dataDims)
 * @param basisFunction The RBF kernel
 * @param inputIDs     IDs of input vertices in inMesh
 * @param inMesh       Input mesh with RBF source vertices
 */
template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
void RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::addWriteDataToCache(const mesh::Vertex &v, const Eigen::VectorXd &load, Eigen::MatrixXd &epsilon, Eigen::MatrixXd &Au,
                                                                         const RADIAL_BASIS_FUNCTION_T &basisFunction, const IndexContainer &inputIDs, const mesh::Mesh &inMesh) const
{
  PRECICE_TRACE();

  // RBF part: for each input vertex j, add phi(||v - x_j||) * load to Au[j].
  // Au accumulates A^T * inputData; each query vertex contributes one row of A^T.
  PRECICE_ASSERT(Au.rows() == static_cast<Eigen::Index>(inputIDs.size()));
  PRECICE_ASSERT(Au.cols() == load.size());

  const auto &out = v.rawCoords();
  for (const auto &j : inputIDs | boost::adaptors::indexed()) {
    const auto &in                = inMesh.vertex(j.value()).rawCoords();
    double      squaredDifference = computeSquaredDifference(out, in);
    double      eval              = basisFunction.evaluate(std::sqrt(squaredDifference));
    Au.row(j.index()) += eval * load.transpose();
  }

  // Polynomial part (SEPARATE mode only):
  // epsilon accumulates V^T * inputData  (polynomial moments of the write data).
  // row 0: constant term contribution  = 1 * load
  // row k: linear term in dimension d  = v.coord[d] * load
  if (epsilon.size() > 0) {
    PRECICE_ASSERT(epsilon.rows() == getNumberOfPolynomials());
    PRECICE_ASSERT(epsilon.cols() == load.size());
    // constant polynomial
    epsilon.row(0) += 1 * load.transpose();
    // the linear contributions — one per active dimension
    int k = 1;
    for (std::size_t d = 0; d < _localActiveAxis.size(); ++d) {
      if (_localActiveAxis[d]) {
        epsilon.row(k) += out[d] * load.transpose();
        ++k;
      }
    }
  }
}

/**
 * EVALUATE CONSERVATIVE CACHE (JIT finalization — conservative)
 *
 * Completes the JIT conservative mapping by solving the linear system using
 * the accumulated buffer `Au` (= A^T * inputData) and the polynomial buffer `epsilon`.
 *
 * This is the single-call finalization step after all per-query `addWriteDataToCache`
 * calls for one time step.
 *
 * Algorithm:
 *   1. out = C^{-1} * Au        (RBF conservative coefficients, "mu" in PETSc notation)
 *   2. (SEPARATE polynomial only):
 *      tau = epsilon - Q^T * out   (polynomial residual)
 *      out -= C_poly^{-T} * tau   (polynomial correction)
 *
 * @param[in,out] epsilon  Polynomial moment buffer (accumulated V^T * writes); overwritten during solve.
 * @param[in]     Au       RBF accumulation buffer (= sum of A[col_j] * load_j)
 * @param[out]    out      Final mapped values at input mesh vertices
 */
template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::evaluateConservativeCache(Eigen::MatrixXd &epsilon, const Eigen::MatrixXd &Au, Eigen::MatrixXd &out) const
{
  // Step 1: Solve for the conservative RBF coefficients.
  // out = C^{-1} * Au  (equivalent to solveConservative's "mu" step)
  out = _decMatrixC.solve(Au);

  if (epsilon.size() > 0) {
    // Step 2a: Compute polynomial residual tau = epsilon - Q^T * mu
    epsilon -= _matrixQ.transpose() * out;
    // Step 2b: Apply polynomial correction: out -= C_poly^{-T} * (-tau) = out + solve(tau)
    out -= static_cast<Eigen::MatrixXd>(_qrMatrixQ.transpose().solve(-epsilon));
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::clear()
{
  _matrixA    = Eigen::MatrixXd();
  _decMatrixC = DecompositionType();
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::Index RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getInputSize() const
{
  return _decMatrixC.cols();
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::Index RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getOutputSize() const
{
  return _matrixA.rows();
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::Index RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getNumberOfPolynomials() const
{
  return 1 + std::count(_localActiveAxis.begin(), _localActiveAxis.end(), true);
}

} // namespace precice::mapping
