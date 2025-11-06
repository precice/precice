#pragma once

#include <mapping/MathHelper.hpp>
#include <mapping/impl/BasisFunctions.hpp>
#include <mesh/Mesh.hpp>

namespace precice::mapping {

struct Sample {
  double radius;
  double error;

  Sample(double radius, double error)
      : radius(radius), error(error)
  {
  }

  Sample()
      : radius(std::numeric_limits<double>::quiet_NaN()), error(std::numeric_limits<double>::max())
  {
  }
};

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
  std::vector<Sample> _samples;

  double _lowerBound;
  double _upperBound;
  bool   _lastSampleWasOptimum;

  Sample _currentOptimum;

  static double estimateMeshResolution(mesh::Mesh &inputMesh);

public:
  ~RBFParameterTuner() = default;
  explicit RBFParameterTuner(mesh::Mesh &inputMesh);

  template <typename IndexContainer, typename Solver>
  Sample optimize(Solver &solver, const IndexContainer &inputIds, const Eigen::VectorXd &inputData);

  bool lastSampleWasOptimum() const;

  double getLastOptimizedRadius() const;
  double getLastOptimizationError() const;

private:
  static bool shouldContinue(const Sample &lowerBound, const Sample &upperBound);
};

inline RBFParameterTuner::RBFParameterTuner(mesh::Mesh &inputMesh)
{
  _lowerBound = estimateMeshResolution(inputMesh);
  _upperBound = std::numeric_limits<double>::quiet_NaN();

  _lastSampleWasOptimum = false;
}

template <typename IndexContainer, typename Solver>
Sample RBFParameterTuner::optimize(Solver &solver, const IndexContainer &inputIds, const Eigen::VectorXd &inputData)
{
  constexpr bool radiusRBF = Solver::BASIS_FUNCTION_T::hasCompactSupport();
  static_assert(radiusRBF, "RBF is not supported by this optimizer, as it does not accept a support-radius."
                            "Currently supported: Compactly supported RBFs and Gaussians.");

  constexpr double increaseSize = 10;
  double           sampleRadius = _lowerBound;

  _lastSampleWasOptimum = false;

  PRECICE_INFO("Starting optimization with lower bound = {:.4e}, upper bound = {:.4e}", _lowerBound, _upperBound);

  solver.rebuildKernelDecomposition(inputIds, sampleRadius);
  auto [error, rcond] = solver.computeErrorEstimate(inputData, inputIds);
  Sample lowerBound   = {sampleRadius, sampleRadius};

  // error is numerically insignificant
  if (error < 1e-15) {
    _lastSampleWasOptimum = true;
    return lowerBound;
  }

  // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
  while (std::isnan(_upperBound)) {
    sampleRadius *= increaseSize;

    solver.rebuildKernelDecomposition(inputIds, sampleRadius);
    auto [error, rcond] = solver.computeErrorEstimate(inputData, inputIds);

    PRECICE_INFO("RBF tuner sample: rad={:.4e}, err={:.4e}, 1/cond={:.4e}", sampleRadius, error, rcond);

    if (rcond < RCOND_TOLERANCE || error >= std::numeric_limits<double>::max()) {
      PRECICE_CHECK(sampleRadius != _lowerBound, "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
      _upperBound = sampleRadius;
      break;
    }
    lowerBound = {sampleRadius, error};
  }
  Sample upperBound = {_upperBound, std::numeric_limits<double>::max()};

  int    i = 0;
  Sample centerSample;

  // bisection algrorithm
  while (shouldContinue(lowerBound, upperBound) && i < MAX_BISEC_ITERATIONS) {
    centerSample.radius = (lowerBound.radius + upperBound.radius) / 2;

    solver.rebuildKernelDecomposition(inputIds, centerSample.radius);
    centerSample.error = std::get<0>(solver.computeErrorEstimate(inputData, inputIds));

    PRECICE_DEBUG("Current interval: [({:.2e},{:.2e}), ({:.2e},{:.2e})], Sample: rad={:.2e}, err={:.2e}",
                 lowerBound.radius, lowerBound.error, upperBound.radius, upperBound.error, centerSample.radius, centerSample.error);

    if (lowerBound.error < upperBound.error || (centerSample.error >= std::numeric_limits<double>::max())) {
      upperBound = centerSample;
    } else {
      lowerBound = centerSample;
    }
    i++;
  }
  _currentOptimum = centerSample.error >= std::numeric_limits<double>::max() ? lowerBound : centerSample;
  _lastSampleWasOptimum = centerSample.radius == _currentOptimum.radius;

  PRECICE_DEBUG("Best sample: rad={:.4e}, err={:.4e}", _currentOptimum.radius, _currentOptimum.error);

  return _currentOptimum;
}

inline double RBFParameterTuner::estimateMeshResolution(mesh::Mesh &inputMesh)
{
  size_t sampleSize = std::min((size_t)3, inputMesh.vertices().size());

  const auto         i0 = inputMesh.vertices().size() / 2;
  const mesh::Vertex x0 = inputMesh.vertices().at(i0);

  const std::vector<VertexID> matches = inputMesh.index().getClosestVertices(x0.getCoords(), sampleSize);

  double h = 0;
  for (int i = 0; i < sampleSize; i++) {
    const mesh::Vertex xi = inputMesh.vertices().at(matches.at(i));
    h += std::sqrt(utils::computeSquaredDifference(xi.rawCoords(), x0.rawCoords()));
  }
  return h / sampleSize;
}

inline bool RBFParameterTuner::shouldContinue(const Sample &lowerBound, const Sample &upperBound)
{
  bool shouldContinue = upperBound.error == std::numeric_limits<double>::max();
  shouldContinue |= upperBound.radius > POS_TOLERANCE * lowerBound.radius;
  shouldContinue |= std::abs(upperBound.error - lowerBound.error) < ERR_TOLERANCE * std::min(upperBound.error, lowerBound.error);
  return shouldContinue;
}

inline bool RBFParameterTuner::lastSampleWasOptimum() const
{
  return _lastSampleWasOptimum;
}

inline double RBFParameterTuner::getLastOptimizedRadius() const {
  return _currentOptimum.radius;
}

inline double RBFParameterTuner::getLastOptimizationError() const {
  return _currentOptimum.error;
}

} // namespace precice::mapping
