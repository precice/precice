#pragma once

#include <mapping/MathHelper.hpp>
#include <mapping/impl/BasisFunctions.hpp>
#include <mesh/Mesh.hpp>

namespace precice::mapping {

struct Sample {
  double pos;
  double error;

  Sample(double pos, double error)
      : pos(pos), error(error)
  {
  }

  Sample()
      : pos(std::numeric_limits<double>::quiet_NaN()), error(std::numeric_limits<double>::infinity())
  {
  }
};

class RBFParameterTuner {

  mutable logging::Logger _log{"mapping::RBFParameterTuner"};

protected:
  Eigen::Index _inSize;

  std::vector<Sample> _samples;

  double _lowerBound;
  double _upperBound;
  bool   _lastSampleWasOptimum;

  static double estimateMeshResolution(mesh::Mesh &inputMesh);

public:
  virtual ~RBFParameterTuner() = default;
  explicit RBFParameterTuner(mesh::Mesh &inputMesh);

  template <typename IndexContainer, typename Solver>
  std::tuple<double, double> optimize(const Solver &solver, const IndexContainer &inputIds, const Eigen::VectorXd &inputData);

  bool lastSampleWasOptimum() const;

private:
  template <typename IndexContainer, typename Solver>
  Sample      optimizeBisection(const Solver &solver, const Eigen::VectorXd &inputData, const IndexContainer &inputIds, double posTolerance, double errorTolerance, int maxIterations);
  static bool shouldContinue(const Sample &lowerBound, const Sample &upperBound, double posTolerance, double errorTolerance);
};

inline RBFParameterTuner::RBFParameterTuner(mesh::Mesh &inputMesh)
{
  _lowerBound = estimateMeshResolution(inputMesh);
  _inSize     = inputMesh.nVertices();
  _upperBound = std::numeric_limits<double>::quiet_NaN();

  _lastSampleWasOptimum = false;
}

template <typename IndexContainer, typename Solver>
std::tuple<double, double> RBFParameterTuner::optimize(const Solver &solver, const IndexContainer &inputIds, const Eigen::VectorXd &inputData)
{
  constexpr bool radiusRBF = Solver::BASIS_FUNCTION_T::hasCompactSupport();
  static_assert(radiusRBF, "RBF is not supported by this optimizer, as it does not accept a support-radius."
                            "Currently supported: Compactly supported RBFs and Gaussians.");

  constexpr double POS_TOLERANCE        = 1.5; // Factor by which the radius is allowed to change before stopping
  constexpr double ERR_TOLERANCE        = 0.5; // Factor by which the error is allowed to change before stopping
  constexpr int    MAX_BISEC_ITERATIONS = 6;   // Number of iterations during the "bisection" step. After finding initial samples

  const Sample bestSample = optimizeBisection(solver, inputData, inputIds, POS_TOLERANCE, ERR_TOLERANCE, MAX_BISEC_ITERATIONS);
  PRECICE_INFO("Best sample: rad={:.4e}, err={:.4e}", bestSample.pos, bestSample.error);

  return {bestSample.pos, bestSample.error};
}

template <typename IndexContainer, typename Solver>
Sample RBFParameterTuner::optimizeBisection(const Solver &solver, const Eigen::VectorXd &inputData, const IndexContainer &inputIds, double posTolerance, double errorTolerance, int maxIterations)
{
  constexpr double increaseSize = 10;
  double           sampleRadius = this->_lowerBound;

  this->_lastSampleWasOptimum = false;

  PRECICE_INFO("Start optimization with lower bound = {:.4e}, upper bound = {:.4e}", this->_lowerBound, this->_upperBound);

  auto [error, rcond] = solver.computeErrorEstimate(inputData, inputIds, sampleRadius);
  Sample lowerBound   = {sampleRadius, sampleRadius};

  // Error is numerically insignificant
  if (error < 1e-15) {
    this->_lastSampleWasOptimum = true;
    return lowerBound;
  }

  // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
  while (std::isnan(this->_upperBound)) {
    sampleRadius *= increaseSize;

    auto [error, rcond] = solver.computeErrorEstimate(inputData, inputIds, sampleRadius);

    PRECICE_INFO("RBF tuner sample: rad={:.4e}, err={:.4e}, 1/cond={:.4e}", sampleRadius, error, rcond);

    if (rcond < 1e-12 || std::isinf(error)) {
      PRECICE_CHECK(sampleRadius != this->_lowerBound, "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
      this->_upperBound = sampleRadius;
      break;
    }
    lowerBound = {sampleRadius, error};
  }
  Sample upperBound = {this->_upperBound, std::numeric_limits<double>::infinity()};

  int    i = 0;
  Sample centerSample;

  while (shouldContinue(lowerBound, upperBound, posTolerance, errorTolerance) && i < maxIterations) {
    centerSample.pos   = (lowerBound.pos + upperBound.pos) / 2;
    centerSample.error = std::get<0>(solver.computeErrorEstimate(inputData, inputIds, centerSample.pos));

    PRECICE_INFO("Current interval: [({:.2e},{:.2e}), ({:.2e},{:.2e})], Sample: rad={:.2e}, err={:.2e}",
                 lowerBound.pos, lowerBound.error, upperBound.pos, upperBound.error, centerSample.pos, centerSample.error);

    if (lowerBound.error < upperBound.error || std::isinf(centerSample.error)) {
      upperBound = centerSample;
    } else {
      lowerBound = centerSample;
    }
    i++;
  }
  const Sample bestSample     = std::isinf(centerSample.error) ? lowerBound : centerSample;
  this->_lastSampleWasOptimum = centerSample.pos == bestSample.pos;

  return bestSample;
}

inline double RBFParameterTuner::estimateMeshResolution(mesh::Mesh &inputMesh)
{
  constexpr int sampleSize = 3;

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

inline bool RBFParameterTuner::shouldContinue(const Sample &lowerBound, const Sample &upperBound, double posTolerance, double errorTolerance)
{
  bool shouldContinue = std::isinf(upperBound.error);
  shouldContinue |= upperBound.pos > posTolerance * lowerBound.pos;
  shouldContinue |= std::abs(upperBound.error - lowerBound.error) < errorTolerance * std::min(upperBound.error, lowerBound.error);
  return shouldContinue;
}

inline bool RBFParameterTuner::lastSampleWasOptimum() const
{
  return _lastSampleWasOptimum;
}

} // namespace precice::mapping
