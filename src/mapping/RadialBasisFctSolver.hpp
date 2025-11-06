#pragma once

#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <boost/range/adaptor/indexed.hpp>
#include <functional>
#include <type_traits>
#include "impl/BasisFunctions.hpp"
#include "impl/BisectionRBFTuner.hpp"
#include "mapping/MathHelper.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/config/MappingConfigurationTypes.hpp"
#include "mapping/impl/BasisFunctions.hpp"
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
  using DecompositionType = std::conditional_t<RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite(), Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>>, Eigen::ColPivHouseholderQR<Eigen::Ref<Eigen::MatrixXd>>>;
  using BASIS_FUNCTION_T  = RADIAL_BASIS_FUNCTION_T;

  /// Deleted since we use Eigen::Ref<Eigen::MatrixXd> for the decomposition type
  RadialBasisFctSolver() = delete;
  /// Deleted since we use Eigen::Ref<Eigen::MatrixXd> for the decomposition type
  RadialBasisFctSolver<BASIS_FUNCTION_T> &operator=(const RadialBasisFctSolver<BASIS_FUNCTION_T>&) = delete;

  /**
   * assembles the system matrices and computes the decomposition of the interpolation matrix
   * inputMesh refers to the mesh where the interpolants are built on, i.e., the input mesh
   * for consistent mappings and the output mesh for conservative mappings
   * outputMesh refers to the mesh where we evaluate the interpolants, i.e., the output mesh
   * consistent mappings and the input mesh for conservative mappings
   */
  template <typename IndexContainer>
  RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, mesh::PtrMesh inputMesh, const IndexContainer &inputIDs,
                       mesh::PtrMesh outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial, MappingConfiguration::AutotuningParams rbfConfig);

  template <typename IndexContainer>
  RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, mesh::PtrMesh inputMesh, const IndexContainer &inputIDs,
                       mesh::PtrMesh outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial); // TODO: refactor and adapt for PUM

  /// Maps the given input data
  template <typename IndexContainer>
  Eigen::VectorXd solveConsistent(Eigen::VectorXd &inputData, const IndexContainer &inputIDs, Polynomial polynomial) const;

  void computeCacheData(Eigen::MatrixXd &inputData, Polynomial polynomial, Eigen::MatrixXd &polyOut, Eigen::MatrixXd &coeffsOut) const;

  /// Maps the given input data
  template <typename IndexContainer>
  Eigen::VectorXd solveConservative(const Eigen::VectorXd &inputData, const IndexContainer &inputIDs, Polynomial polynomial) const;

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

  template <typename IndexConatiner>
  std::tuple<double, double> computeErrorEstimate(const Eigen::VectorXd &inputData, const IndexConatiner &inputIds, double parameter) const;

private:
  mutable precice::logging::Logger _log{"mapping::RadialBasisFctSolver"};

  template <typename IndexConatiner>
  void rebuildKernelDecomposition(const IndexConatiner &inputIds, double parameter) const;

  double evaluateRippaLOOCVerror(const Eigen::VectorXd &lambda) const;

  /// Matrix used for storing the kernel coefficients and the modified coefficients used by the decomposition.
  mutable Eigen::MatrixXd   _decMatrixCoefficients;
  /// Eigen in-place decomposition type of the interpolation matrix _decMatrixCoefficients
  mutable DecompositionType _decMatrixC;

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

  MappingConfiguration::AutotuningParams     _tuningConfig;
  mutable std::unique_ptr<RBFParameterTuner> _tuner;

  mutable int    _iterSinceOptimization = 1;
  mutable double _lastTuningError       = 0;
  mutable double _optimizedRadius       = 0;

  Polynomial _polynomial;

  std::array<bool, 3> _activeAxis;
  mesh::PtrMesh       _inputMesh;
};

// ------- Non-Member Functions ---------

/// given the active axis, computes sets the axis with the lowest spatial expansion to dead
template <typename IndexContainer>
constexpr void reduceActiveAxis(const mesh::PtrMesh &mesh, const IndexContainer &IDs, std::array<bool, 3> &axis)
{
  // make a pair of the axis and the difference
  std::array<std::pair<int, double>, 3> differences;

  // Compute the difference magnitude per direction
  for (std::size_t d = 0; d < axis.size(); ++d) {
    // Ignore dead axis here, i.e., apply the max value such that they are sorted on the last position(s)
    if (axis[d] == false) {
      differences[d] = std::make_pair<int, double>(d, std::numeric_limits<double>::max());
    } else {
      auto res = std::minmax_element(IDs.begin(), IDs.end(), [&](const auto &a, const auto &b) { return mesh->vertex(a).coord(d) < mesh->vertex(b).coord(d); });
      // Check if we are above or below the threshold
      differences[d] = std::make_pair<int, double>(d, std::abs(mesh->vertex(*res.second).coord(d) - mesh->vertex(*res.first).coord(d)));
    }
  }

  std::sort(differences.begin(), differences.end(), [](const auto &d1, const auto &d2) { return d1.second < d2.second; });
  // Disable the axis having the smallest expansion
  axis[differences[0].first] = false;
}

// Fill in the polynomial entries
template <typename IndexContainer>
inline void fillPolynomialEntries(Eigen::MatrixXd &matrix, const mesh::PtrMesh &mesh, const IndexContainer &IDs, Eigen::Index startIndex, std::array<bool, 3> activeAxis)
{
  // Loop over all vertices in the mesh
  for (const auto &i : IDs | boost::adaptors::indexed()) {

    // 1. the constant contribution
    matrix(i.index(), startIndex) = 1.0;

    // 2. the linear contribution
    const auto  &u = mesh->vertex(i.value()).rawCoords();
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
inline void fillKernelEntries(Eigen::MatrixXd &matrixCLU, RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::PtrMesh &inputMesh, const IndexContainer &inputIDs,
                              std::array<bool, 3> activeAxis)
{
  const size_t n = inputIDs.size();
  // Compute RBF matrix entries
  auto         i_iter  = inputIDs.begin();
  Eigen::Index i_index = 0;
  for (; i_iter != inputIDs.end(); ++i_iter, ++i_index) {
    const auto &u       = inputMesh->vertex(*i_iter).rawCoords();
    auto        j_iter  = i_iter;
    auto        j_index = i_index;
    for (; j_iter != inputIDs.end(); ++j_iter, ++j_index) {
      const auto  &v                 = inputMesh->vertex(*j_iter).rawCoords();
      const double squaredDifference = utils::computeSquaredDifference(u, v, activeAxis);
      matrixCLU(i_index, j_index)    = basisFunction.evaluate(std::sqrt(squaredDifference));
    }
  }
  matrixCLU.block(0, 0, n, n).triangularView<Eigen::Lower>() = matrixCLU.block(0, 0, n, n).transpose();
}

template <typename IndexContainer>
inline Eigen::MatrixXd allocateMatrixCLU(const mesh::PtrMesh &inputMesh, const IndexContainer &inputIDs, std::vector<bool> &deadAxis, Polynomial polynomial)
{
  std::array<bool, 3> activeAxis = {false, false, false};
  std::transform(deadAxis.begin(), deadAxis.end(), activeAxis.begin(), [](const auto ax) { return !ax; });

  // Treat the 2D case as 3D case with dead axis
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  const unsigned int dimensions     = 3;
  const unsigned int polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  // Add linear polynom degrees if polynomial requires this
  const auto inputSize = inputIDs.size();
  const auto n         = inputSize + polyparams;

  PRECICE_ASSERT((inputMesh->getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  return Eigen::MatrixXd(n, n);
}

template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
void fillMatrixCLU(Eigen::MatrixXd &matrixCLU, RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::PtrMesh &inputMesh, const IndexContainer &inputIDs,
                                      std::array<bool, 3> activeAxis, Polynomial polynomial)
{
  // Required to fill the poly -> poly entries in the matrix, which remain otherwise untouched
  if (polynomial == Polynomial::ON) {
    matrixCLU.setZero();
  }

  fillKernelEntries(matrixCLU, basisFunction, inputMesh, inputIDs, activeAxis);

  const auto inSize = inputIDs.size();
  const auto polyparams = matrixCLU.cols() - inSize;

  // Add potentially the polynomial contribution in the matrix
  if (polynomial == Polynomial::ON) {
    fillPolynomialEntries(matrixCLU, inputMesh, inputIDs, inSize, activeAxis);
    matrixCLU.block(inSize, 0, polyparams, inSize)                                                  = matrixCLU.transpose().block(inSize, 0, polyparams, inSize);
    matrixCLU.block(inSize, inSize, polyparams, polyparams).template triangularView<Eigen::Lower>() = matrixCLU.block(inSize, inSize, polyparams, polyparams).transpose();
  }
}


template <typename RADIAL_BASIS_FUNCTION_T, typename IndexContainer>
inline Eigen::MatrixXd buildMatrixA(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::PtrMesh &inputMesh, const IndexContainer &inputIDs,
                                    const mesh::PtrMesh &outputMesh, const IndexContainer outputIDs, std::array<bool, 3> activeAxis, Polynomial polynomial)
{
  // Treat the 2D case as 3D case with dead axis
  const unsigned int deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  const unsigned int dimensions     = 3;
  const unsigned int polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  const auto inputSize  = inputIDs.size();
  const auto outputSize = outputIDs.size();
  const auto n          = inputSize + polyparams;

  PRECICE_ASSERT((inputMesh->getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  Eigen::MatrixXd matrixA(outputSize, n);

  // Compute RBF values for matrix A
  for (const auto &i : outputIDs | boost::adaptors::indexed()) {
    const auto &u = outputMesh->vertex(i.value()).rawCoords();
    for (const auto &j : inputIDs | boost::adaptors::indexed()) {
      const auto &v                 = inputMesh->vertex(j.value()).rawCoords();
      double      squaredDifference = utils::computeSquaredDifference(u, v, activeAxis);
      matrixA(i.index(), j.index()) = basisFunction.evaluate(std::sqrt(squaredDifference));
    }
  }

  // Add potentially the polynomial contribution in the matrix
  if (polynomial == Polynomial::ON) {
    fillPolynomialEntries(matrixA, outputMesh, outputIDs, inputSize, activeAxis);
  }
  return matrixA;
}

// ------- Member Functions ---------

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
RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, mesh::PtrMesh inputMesh, const IndexContainer &inputIDs,
                                                                    mesh::PtrMesh outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial)
    : RadialBasisFctSolver(basisFunction, inputMesh, inputIDs, outputMesh, outputIDs, deadAxis, polynomial, {})
{
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, mesh::PtrMesh inputMesh, const IndexContainer &inputIDs,
                                                                    mesh::PtrMesh outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial, MappingConfiguration::AutotuningParams tuningConfig)
    : _decMatrixCoefficients(allocateMatrixCLU(inputMesh, inputIDs, deadAxis, polynomial)), _decMatrixC(_decMatrixCoefficients), _polynomial(polynomial), _activeAxis({false, false, false}), _inputMesh(inputMesh)
{
  PRECICE_ASSERT(!(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite() && polynomial == Polynomial::ON), "The integrated polynomial (polynomial=\"on\") is not supported for the selected radial-basis function. Please select another radial-basis function or change the polynomial configuration.");
  // Convert dead axis vector into an active axis array so that we can handle the reduction more easily
  std::transform(deadAxis.begin(), deadAxis.end(), _activeAxis.begin(), [](const auto ax) { return !ax; });

  _tuningConfig = tuningConfig;

  std::cout << std::endl;
  std::cout << "[ RadialBasisFctSolver() ] _decMatrixC: " << _decMatrixC.rows() << "x" << _decMatrixC.cols() << ", in/outIDs: " << inputIDs.size() << ", " << outputIDs.size() << ", polynomial: " << (polynomial == Polynomial::ON) << std::endl;
  std::cout << std::endl;

  // First, assemble the interpolation matrix and check the invertability
  fillMatrixCLU(_decMatrixCoefficients, basisFunction, inputMesh, inputIDs, _activeAxis, polynomial);
  _decMatrixC.compute(_decMatrixCoefficients);

  if (_tuningConfig.autotuneShape && RADIAL_BASIS_FUNCTION_T::hasCompactSupport()) {
    _tuner = std::make_unique<RBFParameterTuner>(*inputMesh.get());
  } else {
    fillMatrixCLU(_decMatrixCoefficients, basisFunction, inputMesh, inputIDs, _activeAxis, polynomial);
    _decMatrixC.compute(_decMatrixCoefficients);
  }

  PRECICE_CHECK(_tuningConfig.autotuneShape || _decMatrixC.info() == Eigen::ComputationInfo::Success,
                "The interpolation matrix of the RBF mapping from mesh \"{}\" to mesh \"{}\" is not invertable. "
                "This means that the mapping problem is not well-posed. "
                "Please check if your coupling meshes are correct (e.g. no vertices are duplicated) or reconfigure "
                "your basis-function (e.g. reduce the support-radius).",
                inputMesh->getName(), outputMesh->getName());

  // For polynomial on, the algorithm might fail in determining the size of the system
  if (polynomial != Polynomial::ON && computeCrossValidation) {
    // TODO: Disable synchronization
    precice::profiling::Event e("map.rbf.computeLOOCV");
    _inverseDiagonal = utils::computeInverseDiagonal(_decMatrixC);
  }
  // Second, assemble evaluation matrix
  if (_tuningConfig.autotuneShape) {
    _matrixA = buildMatrixA(VolumeSplines(), inputMesh, inputIDs, outputMesh, outputIDs, _activeAxis, polynomial);
  } else {
    _matrixA = buildMatrixA(basisFunction, inputMesh, inputIDs, outputMesh, outputIDs, _activeAxis, polynomial);
  }

  // In case we deal with separated polynomials, we need dedicated matrices for the polynomial contribution
  if (polynomial == Polynomial::SEPARATE) {

    // 4 = 1 + dimensions(3) = maximum number of polynomial parameters
    _localActiveAxis        = _activeAxis;
    unsigned int polyParams = getNumberOfPolynomials();

    do {
      // First, build matrix Q and check for the condition number
      _matrixQ.resize(inputIDs.size(), polyParams);
      fillPolynomialEntries(_matrixQ, inputMesh, inputIDs, 0, _localActiveAxis);

      // Compute the condition number
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(_matrixQ);
      PRECICE_ASSERT(svd.singularValues().size() > 0);
      PRECICE_DEBUG("Singular values in polynomial solver: {}", svd.singularValues());
      const double conditionNumber = svd.singularValues()(0) / std::max(svd.singularValues()(svd.singularValues().size() - 1), math::NUMERICAL_ZERO_DIFFERENCE);
      PRECICE_DEBUG("Condition number: {}", conditionNumber);

      // Disable one axis
      if (conditionNumber > 1e5) {
        reduceActiveAxis(inputMesh, inputIDs, _localActiveAxis);
        polyParams = getNumberOfPolynomials();
      } else {
        break;
      }
    } while (true);

    // allocate and fill matrix V for the outputMesh
    _matrixV.resize(outputIDs.size(), polyParams);
    fillPolynomialEntries(_matrixV, outputMesh, outputIDs, 0, _localActiveAxis);

    // 3. compute decomposition
    _qrMatrixQ = _matrixQ.colPivHouseholderQr();
  }

  if (!_tuningConfig.autotuneShape) {
    // Compute the condition number
    profiling::Event e("map.rbf.condition");

    double rcond = utils::approximateReciprocalConditionNumber(_decMatrixC);
    PRECICE_DEBUG("reciprocal condition number < {}", rcond);
    e.addData("100-log-rcond", static_cast<int>(100 * std::log10(rcond)));
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
Eigen::VectorXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConservative(const Eigen::VectorXd &inputData, const IndexContainer &inputIds, Polynomial polynomial) const
{
  PRECICE_ASSERT((_matrixV.size() > 0 && polynomial == Polynomial::SEPARATE) || _matrixV.size() == 0, _matrixV.size());
  // TODO: Avoid temporary allocations
  // Au is equal to the eta in our PETSc implementation
  PRECICE_ASSERT(inputData.size() == _matrixA.rows());
  Eigen::VectorXd Au = _matrixA.transpose() * inputData;
  PRECICE_ASSERT(Au.size() == _matrixA.cols());

  // mu in the PETSc implementation
  Eigen::VectorXd out = _decMatrixC.solve(Au);

  if (polynomial == Polynomial::SEPARATE) {
    Eigen::VectorXd epsilon = _matrixV.transpose() * inputData;
    PRECICE_ASSERT(epsilon.size() == _matrixV.cols());

    // epsilon = Q^T * mu - epsilon (tau in the PETSc impl)
    epsilon -= _matrixQ.transpose() * out;
    PRECICE_ASSERT(epsilon.size() == _matrixQ.cols());

    // out  = out - solveTranspose tau (sigma in the PETSc impl)
    out -= static_cast<Eigen::VectorXd>(_qrMatrixQ.transpose().solve(-epsilon));
  }
  return out;
}

template <typename RADIAL_BASIS_FUNCTION_T>
inline auto applyKernelToMatrix(const Eigen::Ref<const Eigen::MatrixXd> &matrixExpr, double sampleRadius)
{
  static_assert(RADIAL_BASIS_FUNCTION_T::hasCompactSupport(), "Selected RBF does not support a radius.");

  double parameter = sampleRadius;
  if constexpr (RADIAL_BASIS_FUNCTION_T::requiresRadiusToShapeConversion()) {
    parameter = RADIAL_BASIS_FUNCTION_T::transformRadiusToShape(sampleRadius);
  }
  RADIAL_BASIS_FUNCTION_T kernel(parameter);
  return matrixExpr.unaryExpr([kernel](double x) { return kernel.evaluate(x); });
}

// @todo: change the signature to Eigen::MatrixXd and process all components at once, the solve function of eigen can handle that
template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
Eigen::VectorXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(Eigen::VectorXd &inputData, const IndexContainer &inputIDs, Polynomial polynomial) const
{
  PRECICE_ASSERT((_matrixQ.size() > 0 && polynomial == Polynomial::SEPARATE) || _matrixQ.size() == 0);
  Eigen::VectorXd polynomialContribution;
  // Solve polynomial QR and subtract it from the input data
  if (polynomial == Polynomial::SEPARATE) {
    polynomialContribution = _qrMatrixQ.solve(inputData);
    inputData -= (_matrixQ * polynomialContribution);
  }

  if (_tuningConfig.autotuneShape) {
    if constexpr (RADIAL_BASIS_FUNCTION_T::hasCompactSupport()) {
    // @todo: Data dimension not considered
    const bool optimizeThisIteration = _lastTuningError < 1e-15 || (_tuningConfig.iterationInterval == _iterSinceOptimization);

      // @todo: Add error criterion.
      if (optimizeThisIteration) {
        auto [radius, error] = _tuner->optimize(*this, inputIDs, inputData);
        _optimizedRadius     = radius;
        _lastTuningError     = error;
        if (!_tuner->lastSampleWasOptimum()) {
          rebuildKernelDecomposition(inputIDs, _optimizedRadius);
        }
        PRECICE_INFO("Using: radius={}, LOOCV={}", _optimizedRadius, error);
        _iterSinceOptimization = 1;
      } else {
        PRECICE_INFO("Using radius={} again", _optimizedRadius);
        _iterSinceOptimization++;
      }
    } else {
      // Should not happen, since an error will already be thrown in the configuration.
      PRECICE_ERROR("Parameter optimization is not available for RBFs without a \"support-radius\" configuration.");
    }
  }

  // Integrated polynomial (and separated)
  PRECICE_ASSERT(inputData.size() == _matrixA.cols());
  std::cout << std::endl;
  std::cout << "inputData: " << inputData.rows() << "x" << inputData.cols() << ", _decMatrixC: " << _decMatrixC.rows() << "x" << _decMatrixC.cols() << ", _decMatrixCoefficients: " << _decMatrixCoefficients.rows() << "x" << _decMatrixCoefficients.cols() << std::endl;
  std::cout << std::endl;

  Eigen::VectorXd p = _decMatrixC.solve(inputData);

  if (polynomial != Polynomial::ON && computeCrossValidation) {
    precice::profiling::Event e("map.rbf.evaluateLOOCV");
    PRECICE_INFO("Cross validation error (LOOCV): {}", evaluateRippaLOOCVerror(p));
  }
  PRECICE_ASSERT(p.size() == _matrixA.cols());

  Eigen::VectorXd out;
  if (_tuningConfig.autotuneShape) {
    if constexpr (RADIAL_BASIS_FUNCTION_T::hasCompactSupport()) {
      const int outSize    = _matrixA.rows();
      const int inSize     = inputData.size();
      const int polyParams = _matrixA.cols() - inSize;

      // not yet necessary: integrated polynomial is only supported for SPD functions
      out = applyKernelToMatrix<RADIAL_BASIS_FUNCTION_T>(_matrixA.block(0, 0, outSize, inSize), _optimizedRadius) * p.segment(0, inSize);
      out += _matrixA.block(0, inSize, outSize, polyParams) * p.segment(inSize, polyParams);
    } else {
      // Should not happen, since an error will already be thrown in the configuration.
      PRECICE_ERROR("Parameter optimization is not available for RBFs without a \"support-radius\" configuration.");
    }
  } else {
    out = _matrixA * p;
  }

  // Add the polynomial part again for separated polynomial
  if (polynomial == Polynomial::SEPARATE) {
    out += (_matrixV * polynomialContribution);
  }
  return out;
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexConatiner>
void RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::rebuildKernelDecomposition(const IndexConatiner &inputIds, double parameter) const
{
  static_assert(RADIAL_BASIS_FUNCTION_T::hasCompactSupport(), "RBF does not support a radius-initialization and was still used to instantiate an optimizer.");

  if constexpr (RADIAL_BASIS_FUNCTION_T::requiresRadiusToShapeConversion()) {
    parameter = RADIAL_BASIS_FUNCTION_T::transformRadiusToShape(parameter);
  }
  RADIAL_BASIS_FUNCTION_T kernel(parameter);
  fillMatrixCLU(_decMatrixCoefficients, kernel, _inputMesh, inputIds, _activeAxis, _polynomial);
  _decMatrixC.compute(_decMatrixCoefficients);
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexConatiner>
std::tuple<double, double> RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::computeErrorEstimate(const Eigen::VectorXd &inputData, const IndexConatiner &inputIds, double parameter) const
{
  rebuildKernelDecomposition(inputIds, parameter);

  if (_decMatrixC.info() != Eigen::ComputationInfo::Success) {
    return {std::numeric_limits<double>::infinity(), 0};
  }
  double error = utils::computeRippaLOOCVerror(_decMatrixC, inputData);
  double rcond = utils::approximateReciprocalConditionNumber(_decMatrixC);

  return {error, rcond};
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
    double      squaredDifference = utils::computeSquaredDifference(out, in);
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
    double      squaredDifference = utils::computeSquaredDifference(out, in);
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
  return _matrixA.cols();
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
