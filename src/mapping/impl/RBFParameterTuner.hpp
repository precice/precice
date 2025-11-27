#pragma once

#include <limits>
#include <mapping/MathHelper.hpp>
#include <mapping/impl/BasisFunctions.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Utils.hpp>
#include <optional>

namespace precice::mapping {

/// Solver estimate of the mapping error and the reciprocal condition number of the interpolation matrix.
struct RBFErrorEstimate {
public:
  double error;
  double rcond;

  RBFErrorEstimate(double error, double rcond)
      : error(error), rcond(rcond) {}
};

/// RBF parameter tuner sample consisting of a support radius and its corresponding mapping error estimate.
struct Sample {
public:
  double radius;
  double error;

  Sample(double radius, double error)
      : radius(radius), error(error) {}

  Sample()
      : radius(std::numeric_limits<double>::quiet_NaN()), error(std::numeric_limits<double>::max()) {}
};

/**
 * Parameter optimization of support radius RBFs using a bisection algorithm.
 */
class RBFParameterTuner {

  mutable logging::Logger _log{"mapping::RBFParameterTuner"};

  /// Factor by which the radius is allowed to change before stopping
  static constexpr double POS_TOLERANCE = 1.5;
  /// Factor by which the error is allowed to change before stopping
  static constexpr double ERR_TOLERANCE = 0.5;
  /// Which reciprocal condition number should be considered too low for a good support radius
  static constexpr double RCOND_TOLERANCE = 1e-12;
  /// Number of iterations during the "bisection" step. After finding initial samples
  static constexpr int MAX_BISEC_ITERATIONS = 6;

protected:
  /// Defines the lower bound of the optimization range.
  double _lowerBound;
  /// Defines the upper bound of the optimization range.
  std::optional<double> _upperBound;

  bool   _lastSampleWasOptimum = false;
  Sample _currentOptimum;

  std::string _meshName;

public:
  /**
   * @param[in] inputMesh refers to the mesh where the interpolants are built on,
   * i.e., the input mesh for consistent mappings and the output mesh for conservative mappings.
   *
   * Non const reference because of @ref RBFParameterTuner::estimateMeshResolution(mesh::Mesh &).
   */
  explicit RBFParameterTuner(mesh::Mesh &inputMesh);

  RBFParameterTuner() = default;

  ~RBFParameterTuner() = default;

  /**
   * @brief Optimizes the support radius of an RBF using the bisection method.
   *
   * The bisection method is only guaranteed to find a minimum for convex optimiation problems, however,
   * in the case of RBFs it is possible to find a reasonably good radius of the right order of magnitude.
   *
   * @param[in] solver an RBF mapping solver, see @ref precice::mapping::RadialBasisFctSolver. Expected are the following methods:
   *  - void Solver::rebuildKernelDecomposition(const IndexContainer &inputIds, double radius);
   *  - ErrorEstimate Solver::computeErrorEstimate(const Eigen::VectorXd &inputData, const IndexContainer &inputIds);
   * @param[in] inputIds Index container of the vertices on the input mesh.
   * @param[in] inputData Input values for which the optimizer tries to reduce the mapping error by invoking the computeErrorEstimate() method.
   */
  template <typename IndexContainer, typename Solver>
  Sample optimize(Solver &solver, const mesh::PtrMesh inMesh, const IndexContainer &inputIds, const Eigen::VectorXd &inputData);

  /// Returns true if the last sample tested by the optimizer also corresponds to the found optimum.
  bool lastSampleWasOptimum() const;

  double getLastOptimizedRadius() const;
  double getLastOptimizationError() const;

private:
  /// Check if the bisection algorithm should continue based on @ref POS_TOLERANCE, @ref ERR_TOLERANCE and @ref MAX_BISEC_ITERATIONS
  static bool shouldContinue(const Sample &lowerBound, const Sample &upperBound);
};

inline RBFParameterTuner::RBFParameterTuner(mesh::Mesh &inputMesh)
    : _upperBound(std::nullopt), _lowerBound(mesh::estimateMeshResolution(inputMesh)), _meshName(inputMesh.getName()) {}

template <typename IndexContainer, typename Solver>
Sample RBFParameterTuner::optimize(Solver &solver, const mesh::PtrMesh inMesh, const IndexContainer &inputIds, const Eigen::VectorXd &inputData)
{
  constexpr bool radiusRBF = Solver::BASIS_FUNCTION_T::hasCompactSupport();
  static_assert(radiusRBF, "RBF is not supported by this optimizer, as it does not accept a support-radius."
                           "Currently supported: Compactly supported RBFs and Gaussians.");

  // To determine the _upperBound we exponentially increase the support radius by increaseSize starting with _lowerBound
  // and check if the condition number is too high or the decomposition fails.
  // TODO: maybe make dependent on domain extent.
  constexpr double increaseSize = 10;
  double           sampleRadius = _lowerBound;

  _lastSampleWasOptimum = false;

  PRECICE_INFO("Starting optimization with lower bound = {:.4e}, upper bound = {:.4e} on input mesh: {}.",
               _lowerBound, _upperBound.value_or(std::numeric_limits<double>::quiet_NaN()), _meshName);

  solver.rebuildKernelDecomposition(inMesh, inputIds, sampleRadius);
  RBFErrorEstimate estimate   = solver.computeErrorEstimate(inputData, inputIds);
  Sample           lowerBound = {sampleRadius, estimate.error};

  // error is numerically insignificant: Either zero data or fully explained by the polynomial
  if (estimate.error < 1e-15) {
    _lastSampleWasOptimum = true;
    _currentOptimum       = lowerBound;
    return lowerBound;
  }

  // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
  while (!_upperBound.has_value()) {
    sampleRadius *= increaseSize;

    solver.rebuildKernelDecomposition(inMesh, inputIds, sampleRadius);
    RBFErrorEstimate estimate = solver.computeErrorEstimate(inputData, inputIds);

    PRECICE_INFO("RBF tuner sample: rad={:.4e}, err={:.4e}, 1/cond={:.4e}} on input mesh: {}",
                 sampleRadius, estimate.error, estimate.rcond, _meshName);

    if (estimate.rcond < RCOND_TOLERANCE || estimate.error == std::numeric_limits<double>::max()) {
      PRECICE_CHECK(sampleRadius != _lowerBound, "RBF parameter tuning diverged in first iteration on input mesh {} using support-radius={}."
                                                 "Try using a different basis function, or manually select support radius smaller than {}.",
                    _meshName, sampleRadius, sampleRadius);

      _upperBound = sampleRadius;
      break;
    }
    lowerBound = {sampleRadius, estimate.error};
  }
  Sample upperBound{_upperBound.value(), std::numeric_limits<double>::max()};

  int    i = 0;
  Sample centerSample;

  // bisection algrorithm
  while (shouldContinue(lowerBound, upperBound) && i < MAX_BISEC_ITERATIONS) {
    centerSample.radius = (lowerBound.radius + upperBound.radius) / 2;

    solver.rebuildKernelDecomposition(inMesh, inputIds, centerSample.radius);
    centerSample.error = solver.computeErrorEstimate(inputData, inputIds).error;

    PRECICE_DEBUG("Current interval: [({:.2e},{:.2e}), ({:.2e},{:.2e})], Sample: rad={:.2e}, err={:.2e}, on input mesh: {}.",
                  lowerBound.radius, lowerBound.error, upperBound.radius, upperBound.error, centerSample.radius, centerSample.error, _meshName);

    if (lowerBound.error < upperBound.error || (centerSample.error == std::numeric_limits<double>::max())) {
      upperBound = centerSample;
    } else {
      lowerBound = centerSample;
    }
    i++;
  }
  _currentOptimum       = centerSample.error >= std::numeric_limits<double>::max() ? lowerBound : centerSample;
  _lastSampleWasOptimum = centerSample.radius == _currentOptimum.radius;

  PRECICE_DEBUG("Best sample: rad={:.4e}, err={:.4e} on input mesh: {}.", _currentOptimum.radius, _currentOptimum.error, _meshName);

  return _currentOptimum;
}

inline bool RBFParameterTuner::shouldContinue(const Sample &lowerBound, const Sample &upperBound)
{
  bool shouldContinue = upperBound.error == std::numeric_limits<double>::max();
  shouldContinue |= upperBound.radius > POS_TOLERANCE * lowerBound.radius;
  // TODO: instead of a simple difference consider some logarithmic criterion, since the errors on meshes with different resolutions converge to different values of varying orders of magnitude?
  shouldContinue |= std::abs(upperBound.error - lowerBound.error) < ERR_TOLERANCE * std::min(upperBound.error, lowerBound.error);
  return shouldContinue;
}

inline bool RBFParameterTuner::lastSampleWasOptimum() const
{
  return _lastSampleWasOptimum;
}

inline double RBFParameterTuner::getLastOptimizedRadius() const
{
  return _currentOptimum.radius;
}

inline double RBFParameterTuner::getLastOptimizationError() const
{
  return _currentOptimum.error;
}

} // namespace precice::mapping
