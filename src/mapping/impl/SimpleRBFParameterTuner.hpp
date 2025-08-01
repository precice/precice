#pragma once

#include "mapping/impl/RBFParameterTuner.hpp"
#include <mesh/Mesh.hpp>
#include <mapping/MathHelper.hpp>
#include <mapping/config/MappingConfigurationTypes.hpp>

namespace precice::mapping {

template <typename RBF_T>
class RBFParameterTunerSimple : public RBFParameterTuner<RBF_T> {

  using DecompositionType = typename RBFParameterTuner<RBF_T>::DecompositionType;

  std::vector<Sample> _samples;

  mutable logging::Logger _log{"mapping::SimpleRBFParameterTuner"};

public:
  template <typename IndexContainer>
  RBFParameterTunerSimple(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool, 3> &activeAxis);

  double optimize(const Eigen::VectorXd &inputData) override;

  Sample optimizeIterativeIncrease(const Eigen::VectorXd &inputData, double posTolerance, size_t maxSuccessfulSamples);
  Sample optimizeBisection(const Eigen::VectorXd &inputData, double posTolerance, double errorTolerance, int maxIterations);

private:
  static bool shouldContinue(const Sample &lowerBound, const Sample &upperBound, double posTolerance, double errorTolerance);
};

// Implementation:

template <typename RBF_T>
template <typename IndexContainer>
RBFParameterTunerSimple<RBF_T>::RBFParameterTunerSimple(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool, 3> &activeAxis)
    : RBFParameterTuner<RBF_T>(inputMesh, inputIDs, polynomial, activeAxis)
{
  PRECICE_ASSERT(this->rbfSupportsRadius(), "RBF is not supported by this optimizer, as it does not accept a support-radius."
                                            "Currently supported: Compactly supported RBFs and Gaussians.");
}

template <typename RBF_T>
double RBFParameterTunerSimple<RBF_T>::optimize(const Eigen::VectorXd &inputData)
{
  PRECICE_ASSERT(this->_isInitialized);

  constexpr double POS_TOLERANCE        = 1.5; // Factor by which the radius is allowed to change before stopping
  constexpr double ERR_TOLERANCE        = 0.5; // Factor by which the error is allowed to change before stopping
  constexpr int    MAX_BISEC_ITERATIONS = 6;   // Number of iterations during the "bisection" step. After finding initial samples
  // constexpr int MAX_SUCCESSFUL_SAMPLES = 10;

  Sample bestSample = optimizeBisection(inputData, POS_TOLERANCE, ERR_TOLERANCE, MAX_BISEC_ITERATIONS);
  PRECICE_INFO("Best sample: rad={:.4e}, err={:.4e}", bestSample.pos, bestSample.error);

  return bestSample.pos;
}

template <typename RBF_T>
Sample RBFParameterTunerSimple<RBF_T>::optimizeIterativeIncrease(const Eigen::VectorXd &inputData, double posTolerance, size_t maxSuccessfulSamples)
{
  double increaseSize = 10;
  double sampleRadius = this->_lowerBound;

  // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
  while (increaseSize > posTolerance && _samples.size() <= maxSuccessfulSamples) {

    // exponential increase of support-radius sampling position until failure
    while (true) {

      DecompositionType dec = this->buildKernelDecomposition(sampleRadius);
      if (dec.info() != Eigen::ComputationInfo::Success) {
        PRECICE_INFO("RBF tuner sample: rad={:.4e}, matrix decomposition failed", sampleRadius);
        PRECICE_CHECK(!_samples.empty(), "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
        sampleRadius = _samples.at(_samples.size() - 1).pos;
        break;
      }

      const double rcond = utils::approximateReciprocalConditionNumber(dec);
      const double error = utils::computeRippaLOOCVerror(dec, inputData);

      PRECICE_INFO("RBF tuner sample: rad={:.4e}, err={:.4e}, 1/cond={:.4e}", sampleRadius, error, rcond);

      if (rcond < 1e-13 || std::isnan(error)) {
        PRECICE_CHECK(!_samples.empty(), "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
        sampleRadius = _samples.at(_samples.size() - 1).pos;
        break;
      }

      _samples.emplace_back(sampleRadius, error);
      sampleRadius *= increaseSize;
    }
    increaseSize = std::sqrt(increaseSize);
    sampleRadius *= increaseSize;
  }

  size_t minIdx = 0;
  for (size_t i = 0; i < _samples.size(); i++) {
    if (_samples[i].error < _samples[minIdx].error) {
      minIdx = i;
    }
  }
  return _samples[minIdx];
}

template <typename RBF_T>
bool RBFParameterTunerSimple<RBF_T>::shouldContinue(const Sample &lowerBound, const Sample &upperBound, double posTolerance, double errorTolerance)
{
  bool shouldContinue = std::isnan(upperBound.error);
  shouldContinue |= upperBound.pos > posTolerance * lowerBound.pos;
  shouldContinue |= std::abs(upperBound.error - lowerBound.error) < errorTolerance * std::min(upperBound.error, lowerBound.error);
  return shouldContinue;
}

template <typename RBF_T>
Sample RBFParameterTunerSimple<RBF_T>::optimizeBisection(const Eigen::VectorXd &inputData, double posTolerance, double errorTolerance, int maxIterations)
{
  PRECICE_ASSERT(this->_isInitialized);

  constexpr double increaseSize = 10;
  double           sampleRadius = this->_lowerBound;

  PRECICE_INFO("Start optimization with lower bound = {:.4e}", this->_lowerBound);

  Sample lowerBound = {sampleRadius, std::numeric_limits<double>::quiet_NaN()};
  Sample upperBound = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

  // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
  while (true) {
    sampleRadius *= increaseSize;

    DecompositionType dec = this->buildKernelDecomposition(sampleRadius);
    if (dec.info() != Eigen::ComputationInfo::Success) {
      PRECICE_INFO("RBF tuner sample: rad={:.4e}, matrix decomposition failed", sampleRadius);
      break;
    }

    double rcond = utils::approximateReciprocalConditionNumber(dec);
    double error = utils::computeRippaLOOCVerror(dec, inputData);

    PRECICE_INFO("RBF tuner sample: rad={:.4e}, err={:.4e}, 1/cond={:.4e}", sampleRadius, error, rcond);

    if (rcond < 1e-13 || std::isnan(error)) {
      PRECICE_CHECK(sampleRadius != this->_lowerBound, "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
      break;
    }
    lowerBound = {sampleRadius, error};
  }
  upperBound = {sampleRadius, std::numeric_limits<double>::quiet_NaN()};

  int    i = 0;
  Sample centerSample;

  while (shouldContinue(lowerBound, upperBound, posTolerance, errorTolerance) && i < maxIterations) {
    centerSample.pos   = (lowerBound.pos + upperBound.pos) / 2;
    centerSample.error = utils::computeRippaLOOCVerror(this->buildKernelDecomposition(centerSample.pos), inputData);

    PRECICE_INFO("Current interval: [{:.4e}, {:.4e}], Sample: rad={:.4e}, err={:.4e}", lowerBound.pos, upperBound.pos, centerSample.pos, centerSample.error);

    if (centerSample.error < upperBound.error || std::isnan(upperBound.error) || std::isnan(centerSample.error)) {
      upperBound = centerSample;
    } else {
      lowerBound = centerSample;
    }
    i++;
  }
  return std::isnan(centerSample.error) ? lowerBound : centerSample;
}

} // namespace precice::mapping
