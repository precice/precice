#pragma once

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

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::MatrixXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConservative(const Eigen::MatrixXd &inputData, Polynomial polynomial) const
{
  PRECICE_ASSERT((_matrixV.size() > 0 && polynomial == Polynomial::SEPARATE) || _matrixV.size() == 0, _matrixV.size());
  // TODO: Avoid temporary allocations
  // Au is equal to the eta in our PETSc implementation
  PRECICE_ASSERT(inputData.rows() == _matrixA.rows());
  Eigen::MatrixXd Au = _matrixA.transpose() * inputData;
  PRECICE_ASSERT(Au.rows() == _matrixA.cols());

  // mu in the PETSc implementation
  Eigen::MatrixXd out = _decMatrixC.solve(Au);

  if (polynomial == Polynomial::SEPARATE) {
    Eigen::MatrixXd epsilon = _matrixV.transpose() * inputData;
    PRECICE_ASSERT(epsilon.rows() == _matrixV.cols());

    // epsilon = Q^T * mu - epsilon (tau in the PETSc impl)
    epsilon -= _matrixQ.transpose() * out;
    PRECICE_ASSERT(epsilon.rows() == _matrixQ.cols());

    // out  = out - solveTranspose tau (sigma in the PETSc impl)
    out -= static_cast<Eigen::MatrixXd>(_qrMatrixQ.transpose().solve(-epsilon));
  }
  return out;
}

// @todo: change the signature to Eigen::MatrixXd and process all components at once, the solve function of eigen can handle that
template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::MatrixXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(Eigen::MatrixXd &inputData, Polynomial polynomial) const
{
  PRECICE_ASSERT((_matrixQ.size() > 0 && polynomial == Polynomial::SEPARATE) || _matrixQ.size() == 0);
  Eigen::MatrixXd polynomialContribution;
  // Solve polynomial QR and subtract it from the input data
  if (polynomial == Polynomial::SEPARATE) {
    polynomialContribution = _qrMatrixQ.solve(inputData);
    inputData -= (_matrixQ * polynomialContribution);
  }

  // Integrated polynomial (and separated)
  PRECICE_ASSERT(inputData.rows() == _matrixA.cols());
  Eigen::MatrixXd p = _decMatrixC.solve(inputData);

  if (polynomial != Polynomial::ON && computeCrossValidation) {
    precice::profiling::Event e("map.rbf.evaluateLOOCV");
    PRECICE_INFO("Cross validation error (LOOCV): {}", evaluateRippaLOOCVerror(p));
  }
  PRECICE_ASSERT(p.rows() == _matrixA.cols());
  Eigen::MatrixXd out = _matrixA * p;

  // Add the polynomial part again for separated polynomial
  if (polynomial == Polynomial::SEPARATE) {
    out += (_matrixV * polynomialContribution);
  }
  return out;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::computeCacheData(Eigen::MatrixXd &inputData, Polynomial polynomial, Eigen::MatrixXd &polyOut, Eigen::MatrixXd &coeffsOut) const
{
  PRECICE_ASSERT((_matrixQ.size() > 0 && polynomial == Polynomial::SEPARATE) || _matrixQ.size() == 0);
  // Solve polynomial QR and subtract it from the input data
  if (polynomial == Polynomial::SEPARATE) {
    polyOut = _qrMatrixQ.solve(inputData);
    inputData -= (_matrixQ * polyOut);
  }

  // Integrated polynomial (and separated)
  coeffsOut = _decMatrixC.solve(inputData);
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
Eigen::VectorXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::interpolateAt(const mesh::Vertex &v, const Eigen::MatrixXd &poly, const Eigen::MatrixXd &coeffs,
                                                                             const RADIAL_BASIS_FUNCTION_T &basisFunction, const IndexContainer &inputIDs, const mesh::Mesh &inMesh) const
{
  PRECICE_TRACE();
  // We ignore the dead axis here for the evaluation matrix (matrixA), as the vertex coordinates are zero for a potential 2d case and there is no option in PUM to set them from the user

  // Two cases we have to distinguish here:
  // 1. there is no polynomial given, then the result is out = _matrixA * p;
  Eigen::VectorXd result(coeffs.cols());
  result.setZero();

  // Compute RBF values for matrix A
  const auto &out = v.rawCoords();
  for (const auto &j : inputIDs | boost::adaptors::indexed()) {
    const auto &in                = inMesh.vertex(j.value()).rawCoords();
    double      squaredDifference = computeSquaredDifference(out, in);
    double      eval              = basisFunction.evaluate(std::sqrt(squaredDifference));

    result += eval * coeffs.row(j.index()).transpose();
  }

  // 2. we have a separate polynomial, then we have to add it again here;
  //   out += (_matrixV * polynomialContribution);
  if (poly.size() > 0) {
    // constant polynomial
    result += 1 * poly.row(0).transpose();
    // the linear contributions
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

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
void RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::addWriteDataToCache(const mesh::Vertex &v, const Eigen::VectorXd &load, Eigen::MatrixXd &epsilon, Eigen::MatrixXd &Au,
                                                                        const RADIAL_BASIS_FUNCTION_T &basisFunction, const IndexContainer &inputIDs, const mesh::Mesh &inMesh) const
{
  PRECICE_TRACE();
  // We ignore the dead axis here for the evaluation matrix (matrixA), as the vertex coordinates are zero for a potential 2d case and there is no option in PUM to set them from the user

  // 1. The matrix contribution
  // Compute RBF values for matrix A
  PRECICE_ASSERT(Au.rows() == static_cast<Eigen::Index>(inputIDs.size()));
  PRECICE_ASSERT(Au.cols() == load.size());

  const auto &out = v.rawCoords();
  for (const auto &j : inputIDs | boost::adaptors::indexed()) {
    const auto &in                = inMesh.vertex(j.value()).rawCoords();
    double      squaredDifference = computeSquaredDifference(out, in);
    double      eval              = basisFunction.evaluate(std::sqrt(squaredDifference));

    Au.row(j.index()) += eval * load.transpose();
  }

  // 2. we have a separate polynomial, then we have to add it again here;
  // Eigen::VectorXd epsilon = _matrixV.transpose() * inputData;
  if (epsilon.size() > 0) {
    PRECICE_ASSERT(epsilon.rows() == getNumberOfPolynomials());
    PRECICE_ASSERT(epsilon.cols() == load.size());
    // constant polynomial
    epsilon.row(0) += 1 * load.transpose();
    // the linear contributions
    int k = 1;
    for (std::size_t d = 0; d < _localActiveAxis.size(); ++d) {
      if (_localActiveAxis[d]) {
        epsilon.row(k) += out[d] * load.transpose();
        ++k;
      }
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::evaluateConservativeCache(Eigen::MatrixXd &epsilon, const Eigen::MatrixXd &Au, Eigen::MatrixXd &out) const
{
  // mu in the PETSc implementation
  out = _decMatrixC.solve(Au);

  if (epsilon.size() > 0) {
    epsilon -= _matrixQ.transpose() * out;
    // out  = out - solveTranspose tau (sigma in the PETSc impl)
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
