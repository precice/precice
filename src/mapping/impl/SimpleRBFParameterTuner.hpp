#pragma once

#include "mapping/impl/RBFParameterTuner.hpp"

namespace precice {
namespace mapping {

template <typename RBF_T>
class RBFParameterTunerSimple : public RBFParameterTuner<RBF_T> {

  double              _lowerBound;
  std::vector<Sample> _samples;

  mutable logging::Logger _log{"mapping::RBFParameterTuner"};

public:
  RBFParameterTunerSimple();

  template <typename IndexContainer>
  void   initialize(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool, 3> &activeAxis);
  double optimize(const Eigen::VectorXd &inputData) override;

  Sample optimizeIterativeIncrease(const Eigen::VectorXd &inputData, double posTolerance, size_t maxSuccessfulSamples);
  Sample optimizeBisection(const Eigen::VectorXd &inputData, double posTolerance, int maxIterations);
};

// Implementation:

template <typename RBF_T>
RBFParameterTunerSimple<RBF_T>::RBFParameterTunerSimple()
    : _lowerBound(std::numeric_limits<double>::quiet_NaN())
{
  PRECICE_ASSERT(RBF_T::isStrictlyPositiveDefinite(), "Non SPD RBFs are currently not supported by this optimizer");
  PRECICE_ASSERT(this->rbfSupportsRadius(), "RBF is not supported by this optimizer, as it does not accept a support-radius."
                                            "Currently supported: Compactly supported RBFs, Thin Plate Splines and Gaussians.");
}

template <typename RBF_T>
template <typename IndexContainer>
void RBFParameterTunerSimple<RBF_T>::initialize(const mesh::Mesh &inputMesh, const IndexContainer &inputIDs, const Polynomial &polynomial, const std::array<bool, 3> &activeAxis)
{
  _lowerBound           = this->estimateMeshResolution(inputMesh);
  this->_inSize         = inputIDs.size();
  this->_distanceMatrix = buildMatrixCLU(VolumeSplines(), inputMesh, inputIDs, activeAxis, polynomial);
  this->_isInitialized  = true;
}

template <typename RBF_T>
double RBFParameterTunerSimple<RBF_T>::optimize(const Eigen::VectorXd &inputData)
{
  PRECICE_ASSERT(this->_isInitialized);

  constexpr double POS_TOLERANCE        = 1.5;
  constexpr int    MAX_BISEC_ITERATIONS = 6;
  // constexpr int MAX_SUCCESSFUL_SAMPLES = 10;

  Sample bestSample = optimizeBisection(inputData, POS_TOLERANCE, MAX_BISEC_ITERATIONS);
  PRECICE_INFO("Best sample: rad={:.4e}, err={:.4e}", bestSample.pos, bestSample.error);

  return bestSample.pos;
}

template <typename RBF_T>
Sample RBFParameterTunerSimple<RBF_T>::optimizeIterativeIncrease(const Eigen::VectorXd &inputData, double posTolerance, size_t maxSuccessfulSamples)
{
  double increaseSize = 10;
  double sampleRadius = _lowerBound;

  // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
  while (increaseSize > posTolerance && _samples.size() <= maxSuccessfulSamples) {

    // exponential increase of support-radius sampling position until failure
    while (true) {

      Eigen::LLT<Eigen::MatrixXd> llt = this->buildKernelLLT(sampleRadius);
      if (llt.info() != Eigen::ComputationInfo::Success) {
        PRECICE_INFO("RBF tuner sample: rad={:.4e}, matrix decomposition failed", sampleRadius);
        PRECICE_CHECK(!_samples.empty(), "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
        sampleRadius = _samples.at(_samples.size() - 1).pos;
        break;
      }

      const double rcond = utils::approximateReciprocalConditionNumber(llt);
      const double error = utils::computeRippaLOOCVerror(llt, inputData);

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
Sample RBFParameterTunerSimple<RBF_T>::optimizeBisection(const Eigen::VectorXd &inputData, double posTolerance, int maxIterations)
{
  PRECICE_ASSERT(this->_isInitialized);

  constexpr double increaseSize = 10;
  double           sampleRadius = _lowerBound;

  PRECICE_INFO("Start optimization with lower bound = {:.4e}", _lowerBound);

  Sample lowerBound = {sampleRadius, std::numeric_limits<double>::quiet_NaN()};
  Sample upperBound = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

  // collect samples with exponential growth and decrease growth rate until the support radius is "good enough"
  while (true) {
    sampleRadius *= increaseSize;

    Eigen::LLT<Eigen::MatrixXd> llt = this->buildKernelLLT(sampleRadius);
    if (llt.info() != Eigen::ComputationInfo::Success) {
      PRECICE_INFO("RBF tuner sample: rad={:.4e}, matrix decomposition failed", sampleRadius);
      break;
    }

    double rcond = utils::approximateReciprocalConditionNumber(llt);
    double error = utils::computeRippaLOOCVerror(llt, inputData);

    PRECICE_INFO("RBF tuner sample: rad={:.4e}, err={:.4e}, 1/cond={:.4e}", sampleRadius, error, rcond);

    if (rcond < 1e-13 || std::isnan(error)) {
      PRECICE_CHECK(sampleRadius != _lowerBound, "Parameter tuning failed in first iteration using support-radius={}", sampleRadius);
      break;
    }
    lowerBound = {sampleRadius, error};
  }
  upperBound = {sampleRadius, std::numeric_limits<double>::quiet_NaN()};

  int    i = 0;
  Sample centerSample;

  while ((std::isnan(upperBound.error) || upperBound.pos > posTolerance * lowerBound.pos) && i < maxIterations) {
    centerSample.pos   = (lowerBound.pos + upperBound.pos) / 2;
    centerSample.error = utils::computeRippaLOOCVerror(this->buildKernelLLT(centerSample.pos), inputData);

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

} // namespace mapping
} // namespace precice