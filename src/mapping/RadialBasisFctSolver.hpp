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
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"
#include "impl/RBFParameterTuner.hpp"
#include "impl/RBFMatrixOperations.hpp"

namespace precice {
namespace mapping {

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
                       const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial, MappingConfiguration::RBFOptional rbfConfig);

  template <typename IndexContainer>
  RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                       const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial); // TODO: refactor and adapt for PUM

  /// Maps the given input data
  Eigen::VectorXd solveConsistent(Eigen::VectorXd &inputData, Polynomial polynomial) const;

  /// Maps the given input data
  Eigen::VectorXd solveConservative(const Eigen::VectorXd &inputData, Polynomial polynomial) const;

  // Clear all stored matrices
  void clear();

  // Returns the size of the input data
  Eigen::Index getInputSize() const;

  // Returns the size of the input data
  Eigen::Index getOutputSize() const;

  void setClusterRadius(double r)
  {
    this->clusterRadius = r;
  }
  // BO_PARAM(std::size_t, dim_in, 1);
  // number of dimensions of the result (res.size())
  // BO_PARAM(std::size_t, dim_out, 1);
private:
  mutable precice::logging::Logger _log{"mapping::RadialBasisFctSolver"};

  bool _autotuneShape; // TODO: add documentation

  double evaluateRippaLOOCVerror(const Eigen::VectorXd &lambda) const;
  /// Decomposition of the interpolation matrix
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
  mutable Eigen::MatrixXd _matrixA;

  mutable RBFParameterTunerSimple<RADIAL_BASIS_FUNCTION_T> _tuner;

  // TODO: Won't work with global RBF, as we set the minimum in the SphericalVertexCLuster as the (half) cluster radius or similar
  double clusterRadius = std::numeric_limits<double>::quiet_NaN();

  bool computeCrossValidation = false;
};


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
                                                                      : RadialBasisFctSolver(basisFunction, inputMesh, inputIDs, outputMesh, outputIDs, deadAxis, polynomial, {})
{ }

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::RadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                                                                    const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial, MappingConfiguration::RBFOptional rbfConfig)
{
  PRECICE_ASSERT(!(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite() && polynomial == Polynomial::ON), "The integrated polynomial (polynomial=\"on\") is not supported for the selected radial-basis function. Please select another radial-basis function or change the polynomial configuration.");
  // Convert dead axis vector into an active axis array so that we can handle the reduction more easily
  std::array<bool, 3> activeAxis({{false, false, false}});
  std::transform(deadAxis.begin(), deadAxis.end(), activeAxis.begin(), [](const auto ax) { return !ax; });

  _autotuneShape = rbfConfig.autotuneShape;

  fmt::println("\nAUTOTUNE SHAPE: {}", _autotuneShape);

  // First, assemble the interpolation matrix and check the invertability
  bool decompositionSuccessful = false;
  if constexpr (RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite()) {
    if (_autotuneShape) {
      _tuner.initialize(inputMesh, inputIDs, polynomial, activeAxis);
    } else {
      _decMatrixC             = buildMatrixCLU(basisFunction, inputMesh, inputIDs, activeAxis, polynomial).llt();
      decompositionSuccessful = _decMatrixC.info() == Eigen::ComputationInfo::Success;
    }
  } else { // TODO: Auto tuning not supported
    _decMatrixC             = buildMatrixCLU(basisFunction, inputMesh, inputIDs, activeAxis, polynomial).colPivHouseholderQr();
    decompositionSuccessful = _decMatrixC.isInvertible();
  }

  // PRECICE_CHECK(decompositionSuccessful,
  //               "The interpolation matrix of the RBF mapping from mesh \"{}\" to mesh \"{}\" is not invertable. "
  //               "This means that the mapping problem is not well-posed. "
  //               "Please check if your coupling meshes are correct (e.g. no vertices are duplicated) or reconfigure "
  //               "your basis-function (e.g. reduce the support-radius).",
  //               inputMesh.getName(), outputMesh.getName());

  // For polynomial on, the algorithm might fail in determining the size of the system
  if (polynomial != Polynomial::ON && computeCrossValidation) {
    // TODO: Disable synchronization
    precice::profiling::Event e("map.rbf.computeLOOCV");
    _inverseDiagonal = computeInverseDiagonal(_decMatrixC);
  }
  // Second, assemble evaluation matrix
  if (_autotuneShape) {
    _matrixA = buildMatrixA(VolumeSplines(), inputMesh, inputIDs, outputMesh, outputIDs, activeAxis, polynomial); // TODO: necessary?
  } else {
    _matrixA = buildMatrixA(basisFunction, inputMesh, inputIDs, outputMesh, outputIDs, activeAxis, polynomial);
  }

  // In case we deal with separated polynomials, we need dedicated matrices for the polynomial contribution
  if (polynomial == Polynomial::SEPARATE) {

    // 4 = 1 + dimensions(3) = maximum number of polynomial parameters
    auto         localActiveAxis = activeAxis;
    unsigned int polyParams      = 4 - std::count(localActiveAxis.begin(), localActiveAxis.end(), false);

    do {
      // First, build matrix Q and check for the condition number
      _matrixQ.resize(inputIDs.size(), polyParams);
      fillPolynomialEntries(_matrixQ, inputMesh, inputIDs, 0, localActiveAxis);

      // Compute the condition number
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(_matrixQ);
      PRECICE_ASSERT(svd.singularValues().size() > 0);
      PRECICE_DEBUG("Singular values in polynomial solver: {}", svd.singularValues());
      const double conditionNumber = svd.singularValues()(0) / std::max(svd.singularValues()(svd.singularValues().size() - 1), math::NUMERICAL_ZERO_DIFFERENCE);
      PRECICE_DEBUG("Condition number: {}", conditionNumber);

      // Disable one axis
      if (conditionNumber > 1e5) {
        reduceActiveAxis(inputMesh, inputIDs, localActiveAxis);
        polyParams = 4 - std::count(localActiveAxis.begin(), localActiveAxis.end(), false);
      } else {
        break;
      }
    } while (true);

    // allocate and fill matrix V for the outputMesh
    _matrixV.resize(outputIDs.size(), polyParams);
    fillPolynomialEntries(_matrixV, outputMesh, outputIDs, 0, localActiveAxis);
    // 3. compute decomposition
    _qrMatrixQ = _matrixQ.colPivHouseholderQr();
  }

  if (!_autotuneShape) {
    // Compute the condition number
    profiling::Event e("map.rbf.condition");
    double cond = approximateConditionNumber(_decMatrixC);
    fmt::print("\nCONDITION <= {}\n", cond);
    e.addData("condition", static_cast<int>(cond));
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::VectorXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConservative(const Eigen::VectorXd &inputData, Polynomial polynomial) const
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
Eigen::VectorXd RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(Eigen::VectorXd &inputData, Polynomial polynomial) const
{
  PRECICE_ASSERT((_matrixQ.size() > 0 && polynomial == Polynomial::SEPARATE) || _matrixQ.size() == 0);
  Eigen::VectorXd polynomialContribution;
  // Solve polynomial QR and subtract it from the input data
  if (polynomial == Polynomial::SEPARATE) {
    polynomialContribution = _qrMatrixQ.solve(inputData);
    inputData -= (_matrixQ * polynomialContribution);
  }

  if (_autotuneShape) {
    //std::cout << "_meshResolution=" << _meshResolution << "\n";
    //std::cout << "clusterRadius=" << clusterRadius << "\n";

    // if (_tuner.parameterError() > )
    fmt::println("start");
    double optRadius = _tuner.optimize(inputData);

    if constexpr (RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite()) {
      RADIAL_BASIS_FUNCTION_T kernel(optRadius);
      _decMatrixC = _tuner.getInterpolationMatrix().llt();
      _matrixA    = _matrixA.unaryExpr([&kernel] (double x) { return kernel.evaluate(x); });
    } else {
      PRECICE_ASSERT(false, "Not supported.");
    }
  }

  // Integrated polynomial (and separated)
  PRECICE_ASSERT(inputData.size() == _matrixA.cols());
  Eigen::VectorXd p = _decMatrixC.solve(inputData);

  if (polynomial != Polynomial::ON && computeCrossValidation) {
    precice::profiling::Event e("map.rbf.evaluateLOOCV");
    PRECICE_INFO("Cross validation error (LOOCV): {}", evaluateRippaLOOCVerror(p));
  }
  PRECICE_ASSERT(p.size() == _matrixA.cols());

  Eigen::VectorXd out = _matrixA * p;

  // Add the polynomial part again for separated polynomial
  if (polynomial == Polynomial::SEPARATE) {
    out += (_matrixV * polynomialContribution);
  }
  return out;
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
} // namespace mapping
} // namespace precice
