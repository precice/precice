#pragma once

#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <boost/range/adaptor/indexed.hpp>
#include "impl/BasisFunctions.hpp"
#include "impl/RBFParameterTuner.hpp"
#include "mapping/MathHelper.hpp"
#include "mapping/RadialBasisFctSolver.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/config/MappingConfigurationTypes.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/Mesh.hpp"
#include "profiling/Event.hpp"

namespace precice::mapping {

/**
 * This is a solver similar to \ref RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T> which implements a consistent mapping for which the RBF support radius parameter is automatically tuned.
 *
 * Conservative parameter tuning is not available since most error metrics decrease alongside a decreasing support radius making minimization difficult.
 * For such mappings, conservativeness can be ensured using a polynomial and interpolated values can be distributed on the mesh with a simple Volume Spline.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class AutoTunedRBFSolver : public RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T> {

  mutable precice::logging::Logger _log{"mapping::AutoTunedRBFSolver"};

public:
  using BASIS_FUNCTION_T = RADIAL_BASIS_FUNCTION_T;
  using AutotuningParams = MappingConfiguration::AutotuningParams;

  template <typename IndexContainer>
  Eigen::VectorXd solveConsistent(Eigen::VectorXd &inputData, const IndexContainer &inputIDs, Polynomial polynomial);

  template <typename IndexContainer>
  AutoTunedRBFSolver(RADIAL_BASIS_FUNCTION_T basisFunction, mesh::PtrMesh inputMesh, const IndexContainer &inputIDs,
                     mesh::PtrMesh outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial, AutotuningParams rbfTunerConfig);

  AutoTunedRBFSolver() = default;

  /**
   * Rebuilds the matrix decomposition based for a given support radius by updating \ref _decMatrixC and \ref _decMatrixCoefficients.
   * @param[in] parameter RBF support radius.
   */
  template <typename IndexContainer>
  void rebuildKernelDecomposition(const mesh::PtrMesh inMesh, const IndexContainer &inputIds, double parameter);

  /**
   * Returns an error estimate containing the LOOCV error and an approximation of the reciprocal condition number
   * based on the current decomposition \ref _decMatrixC.
   */
  template <typename IndexContainer>
  RBFErrorEstimate computeErrorEstimate(const Eigen::VectorXd &inputData, const IndexContainer &inputIds) const;

private:
  AutotuningParams  _rbfTunerConfig;
  RBFParameterTuner _tuner;

  /// The tuner needs access to the global input mesh in \ref solveConsistent to recompute the matrix decomposition.
  /// No ownership is taken and the tuner is not available for the gather-scatter parallelism of \ref mapping::RadialBasisFctMapping to avoid copies of the global mesh.
  std::weak_ptr<mesh::Mesh> _inputMesh;

  Polynomial _polynomial;

  int _iterSinceOptimization = 1;
};

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

  // Integrated polynomial (and separated)
  PRECICE_ASSERT(inputData.size() == _matrixA.cols());
  Eigen::VectorXd p = _decMatrixC->solve(inputData);

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
template <typename IndexContainer>
AutoTunedRBFSolver<RADIAL_BASIS_FUNCTION_T>::AutoTunedRBFSolver(RADIAL_BASIS_FUNCTION_T basisFunction, mesh::PtrMesh inputMesh, const IndexContainer &inputIDs,
                                                                mesh::PtrMesh outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis,
                                                                Polynomial polynomial, AutotuningParams rbfTunerConfig)
    : RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>(basisFunction, inputMesh, inputIDs, deadAxis, polynomial),
      _rbfTunerConfig(rbfTunerConfig),
      _inputMesh(inputMesh),
      _tuner(*inputMesh.get()),
      _polynomial(polynomial)
{
  // For polynomial on, the algorithm might fail in determining the size of the system
  if (polynomial != Polynomial::ON && this->computeCrossValidation) {
    // TODO: Disable synchronization
    precice::profiling::Event e("map.rbf.computeLOOCV");
    this->_inverseDiagonal = utils::computeInverseDiagonal(*(this->_decMatrixC));
  }
  // Allocate and initialize the evaluation matrix. We use VolumeSplines, since the actual entries depend on the result of the tuner.
  this->_matrixA = buildMatrixA(VolumeSplines(), inputMesh, inputIDs, outputMesh, outputIDs, this->_activeAxis, polynomial);

  // In case we deal with separated polynomials, we need dedicated matrices for the polynomial contribution
  if (polynomial == Polynomial::SEPARATE) {
    this->configureSeparatePolynomial(inputMesh, inputIDs, outputMesh, outputIDs);
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
void AutoTunedRBFSolver<RADIAL_BASIS_FUNCTION_T>::rebuildKernelDecomposition(const mesh::PtrMesh inMesh, const IndexContainer &inputIds, double parameter)
{
  static_assert(RADIAL_BASIS_FUNCTION_T::hasCompactSupport(), "RBF does not support a radius-initialization and was still used to instantiate an optimizer.");

  if constexpr (RADIAL_BASIS_FUNCTION_T::requiresRadiusToShapeConversion()) {
    parameter = RADIAL_BASIS_FUNCTION_T::transformRadiusToShape(parameter);
  }
  RADIAL_BASIS_FUNCTION_T kernel(parameter);
  fillMatrixCLU(this->_decMatrixCoefficients, kernel, inMesh, inputIds, this->_activeAxis, this->_polynomial);
  this->_decMatrixC->compute(this->_decMatrixCoefficients);
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
RBFErrorEstimate AutoTunedRBFSolver<RADIAL_BASIS_FUNCTION_T>::computeErrorEstimate(const Eigen::VectorXd &inputData, const IndexContainer &inputIds) const
{
  if (this->_decMatrixC->info() != Eigen::ComputationInfo::Success) {
    return {std::numeric_limits<double>::max(), 0};
  }
  double error = utils::computeRippaLOOCVerror(*(this->_decMatrixC), inputData);
  double rcond = utils::approximateReciprocalConditionNumber(*(this->_decMatrixC));

  return {error, rcond};
}

// @todo: change the signature to Eigen::MatrixXd and process all components at once, the solve function of eigen can handle that
template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
Eigen::VectorXd AutoTunedRBFSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(Eigen::VectorXd &inputData, const IndexContainer &inputIDs, Polynomial polynomial)
{
  PRECICE_ASSERT((this->_matrixQ.size() > 0 && polynomial == Polynomial::SEPARATE) || this->_matrixQ.size() == 0);
  Eigen::VectorXd polynomialContribution;
  // Solve polynomial QR and subtract it from the input data
  if (polynomial == Polynomial::SEPARATE) {
    polynomialContribution = this->_qrMatrixQ.solve(inputData);
    inputData -= (this->_matrixQ * polynomialContribution);
  }

  if (_rbfTunerConfig.autotuneShape) {
    if constexpr (RADIAL_BASIS_FUNCTION_T::hasCompactSupport()) {

      const bool isProbablyZeroData = _tuner.getLastOptimizationError() < 1e-15;
      // @todo: Data dimension not considered
      const bool isOptimizationInterval   = _rbfTunerConfig.iterationInterval == _iterSinceOptimization;
      const bool isOptimizedAndShouldOnce = _rbfTunerConfig.optimizeOnce() && _tuner.getLastOptimizationError() != std::numeric_limits<double>::max();

      // @todo: Add error criterion.
      if (isProbablyZeroData || isOptimizationInterval || !isOptimizedAndShouldOnce) {
        Sample sample = _tuner.optimize(*this, this->_inputMesh.lock(), inputIDs, inputData);
        if (!_tuner.lastSampleWasOptimum()) {
          this->rebuildKernelDecomposition(this->_inputMesh.lock(), inputIDs, _tuner.getLastOptimizedRadius());
        }
        PRECICE_INFO("Using: radius={}, LOOCV={}", _tuner.getLastOptimizedRadius(), sample.error);
        _iterSinceOptimization = 1;
      } else {
        PRECICE_INFO("Using radius={} again", _tuner.getLastOptimizedRadius());
        _iterSinceOptimization++;
      }
    } else {
      // Should not happen, since an error will already be thrown in the configuration.
      PRECICE_ERROR("Parameter optimization is not available for RBFs without a \"support-radius\" configuration.");
    }
  }

  // Integrated polynomial (and separated)
  PRECICE_ASSERT(inputData.size() == this->_matrixA.cols());
  Eigen::VectorXd p = this->_decMatrixC->solve(inputData);

  // if (polynomial != Polynomial::ON && computeCrossValidation) {
  //   precice::profiling::Event e("map.rbf.evaluateLOOCV");
  //   PRECICE_INFO("Cross validation error (LOOCV): {}", evaluateRippaLOOCVerror(p));
  // }
  PRECICE_ASSERT(p.size() == this->_matrixA.cols());

  Eigen::VectorXd out;
  if (_rbfTunerConfig.autotuneShape) {
    if constexpr (RADIAL_BASIS_FUNCTION_T::hasCompactSupport()) {
      const Eigen::Index outSize    = this->_matrixA.rows();
      const Eigen::Index inSize     = inputData.size();
      const Eigen::Index polyParams = this->_matrixA.cols() - inSize;

      // not yet necessary: integrated polynomial is only supported for SPD functions
      out = applyKernelToMatrix<RADIAL_BASIS_FUNCTION_T>(this->_matrixA.block(0, 0, outSize, inSize), _tuner.getLastOptimizedRadius()) * p.segment(0, inSize);
      out += this->_matrixA.block(0, inSize, outSize, polyParams) * p.segment(inSize, polyParams);
    } else {
      // Should not happen, since an error will already be thrown in the configuration.
      PRECICE_ERROR("Parameter optimization is not available for RBFs without a \"support-radius\" configuration.");
    }
  } else {
    out = this->_matrixA * p;
  }

  // Add the polynomial part again for separated polynomial
  if (polynomial == Polynomial::SEPARATE) {
    out += (this->_matrixV * polynomialContribution);
  }
  return out;
}

}; // namespace precice::mapping
