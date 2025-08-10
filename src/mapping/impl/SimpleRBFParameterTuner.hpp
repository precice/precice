#pragma once

#include <mapping/MathHelper.hpp>
#include <mapping/config/MappingConfigurationTypes.hpp>
#include <mesh/Mesh.hpp>
#include "mapping/impl/RBFParameterTuner.hpp"

namespace precice::mapping {

template <typename Solver>
class RBFParameterTunerSimple : public RBFParameterTuner<Solver> {

  mutable logging::Logger _log{"mapping::SimpleRBFParameterTuner"};

  std::vector<Sample> _samples;

public:
  explicit RBFParameterTunerSimple(const mesh::Mesh &inputMesh);

  double optimize(const Solver &solver, const Eigen::VectorXd &inputData) override;
  Sample optimizeBisection(const Solver &solver, const Eigen::VectorXd &inputData, double posTolerance, double errorTolerance, int maxIterations);

private:
  static bool shouldContinue(const Sample &lowerBound, const Sample &upperBound, double posTolerance, double errorTolerance);
};

template <typename Solver>
RBFParameterTunerSimple<Solver>::RBFParameterTunerSimple(const mesh::Mesh &inputMesh)
    : RBFParameterTuner<Solver>(inputMesh)
{
}

template <typename Solver>
double RBFParameterTunerSimple<Solver>::optimize(const Solver &solver, const Eigen::VectorXd &inputData)
{
  constexpr double POS_TOLERANCE        = 1.5; // Factor by which the radius is allowed to change before stopping
  constexpr double ERR_TOLERANCE        = 0.5; // Factor by which the error is allowed to change before stopping
  constexpr int    MAX_BISEC_ITERATIONS = 6;   // Number of iterations during the "bisection" step. After finding initial samples

  Sample bestSample = optimizeBisection(solver, inputData, POS_TOLERANCE, ERR_TOLERANCE, MAX_BISEC_ITERATIONS);
  PRECICE_INFO("Best sample: rad={:.4e}, err={:.4e}", bestSample.pos, bestSample.error);

  return bestSample.pos;
}

template <typename RBF_T>
bool RBFParameterTunerSimple<RBF_T>::shouldContinue(const Sample &lowerBound, const Sample &upperBound, double posTolerance, double errorTolerance)
{
  bool shouldContinue = std::isnan(upperBound.error);
  shouldContinue |= upperBound.pos > posTolerance * lowerBound.pos;
  shouldContinue |= std::abs(upperBound.error - lowerBound.error) < errorTolerance * std::min(upperBound.error, lowerBound.error);
  return shouldContinue;
}

template <typename Solver>
Sample RBFParameterTunerSimple<Solver>::optimizeBisection(const Solver &solver, const Eigen::VectorXd &inputData, double posTolerance, double errorTolerance, int maxIterations)
{
  constexpr double increaseSize = 10;
  double           sampleRadius = this->_lowerBound;

  PRECICE_INFO("Start optimization with lower bound = {:.4e}", this->_lowerBound);

  Sample lowerBound = {sampleRadius, std::numeric_limits<double>::quiet_NaN()};

  // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
  while (true) {
    sampleRadius *= increaseSize;

    auto [error, rcond] = solver.computeErrorEstimate(inputData, sampleRadius);

    PRECICE_INFO("RBF tuner sample: rad={:.4e}, err={:.4e}, 1/cond={:.4e}", sampleRadius, error, rcond);

    if (rcond < 1e-13 || std::isinf(error)) {
      PRECICE_CHECK(sampleRadius != this->_lowerBound, "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
      break;
    }
    lowerBound = {sampleRadius, error};
  }
  Sample upperBound = {sampleRadius, std::numeric_limits<double>::infinity()};

  int    i = 0;
  Sample centerSample;

  while (shouldContinue(lowerBound, upperBound, posTolerance, errorTolerance) && i < maxIterations) {
    centerSample.pos   = (lowerBound.pos + upperBound.pos) / 2;
    centerSample.error = std::get<0>(solver.computeErrorEstimate(inputData, centerSample.pos));

    PRECICE_INFO("Current interval: [{:.4e}, {:.4e}], Sample: rad={:.4e}, err={:.4e}", lowerBound.pos, upperBound.pos, centerSample.pos, centerSample.error);

    if (lowerBound.error < upperBound.error || std::isinf(centerSample.error)) {
      upperBound = centerSample;
    } else {
      lowerBound = centerSample;
    }
    i++;
  }
  const Sample optimum        = std::isnan(centerSample.error) ? lowerBound : centerSample;
  this->_lastSampleWasOptimum = centerSample.pos == optimum.pos;

  return optimum;
}

} // namespace precice::mapping
