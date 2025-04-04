#include "acceleration/impl/ResidualSumPreconditioner.hpp"
#include <algorithm>
#include <cmath>
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::acceleration::impl {

ResidualSumPreconditioner::ResidualSumPreconditioner(
    int  maxNonConstTimeWindows,
    bool preconditionerUpdateOnThreshold)
    : Preconditioner(maxNonConstTimeWindows),
      _preconditionerUpdateOnThreshold(preconditionerUpdateOnThreshold)
{
}

void ResidualSumPreconditioner::initialize(std::vector<size_t> &svs)
{
  PRECICE_TRACE();
  Preconditioner::initialize(svs);

  _residualSum.resize(_subVectorSizes.size(), 0.0);
  _previousResidualSum.resize(_subVectorSizes.size(), 0.0);
}

void ResidualSumPreconditioner::_update_(bool                   timeWindowComplete,
                                         const Eigen::VectorXd &oldValues,
                                         const Eigen::VectorXd &res)
{
  if (timeWindowComplete) {
    _firstTimeWindow = false;
    std::fill(_residualSum.begin(), _residualSum.end(), 0.0);
    return;
  }

  std::vector<double> norms(_subVectorSizes.size(), 0.0);

  double sum = 0.0;

  for (size_t k = 0; k < _subVectorSizes.size(); k++) {
    Eigen::VectorXd part = res.segment(_subVectorOffsets[k], _subVectorSizes[k]);
    norms[k]             = utils::IntraComm::dot(part, part);
    sum += norms[k];
    norms[k] = std::sqrt(norms[k]);
  }
  sum = std::sqrt(sum);
  PRECICE_WARN_IF(
      math::equals(sum, 0.0),
      "All residual sub-vectors in the residual-sum preconditioner are numerically zero ( sum = {}). "
      "This indicates that the data values exchanged between two successive iterations did not change. "
      "The simulation may be unstable, e.g. produces NAN values. Please check the data values exchanged "
      "between the solvers is not identical between iterations. The preconditioner scaling factors were "
      "not updated in this iteration and the scaling factors determined in the previous iteration were used.",
      sum);

  for (size_t k = 0; k < _subVectorSizes.size(); k++) {
    if (sum > math::NUMERICAL_ZERO_DIFFERENCE)
      _residualSum[k] += norms[k] / sum;

    PRECICE_WARN_IF(
        math::equals(_residualSum[k], 0.0),
        "A sub-vector in the residual-sum preconditioner became numerically zero ( sub-vector = {}). "
        "If this occurred in the second iteration and the initial-relaxation factor is equal to 1.0, "
        "check if the coupling data values of one solver is zero in the first iteration. "
        "The preconditioner scaling factors were not updated for this iteration and the scaling factors "
        "determined in the previous iteration were used.",
        _residualSum[k]);
  }

  // Determine if weights needs to be reset
  // if _preconditionerUpdateOnThreshold is true, the weights are reset only if the ratio of the new scaling weight to the previous residual sum has changed significantly
  bool resetWeights = _firstTimeWindow || !_preconditionerUpdateOnThreshold;
  if (!resetWeights) {
    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      double resSum = _residualSum[k];
      if (math::equals(resSum, 0.0)) {
        continue; // These will be ignored when resetting the weights
      }
      // Check if the ratio of the new scaling weight to the previous residual sum
      // has changed significantly, either exceeding the threshold of 10.0
      // or dropping below its inverse.
      double factor = _previousResidualSum[k] / resSum;
      if ((factor > 10.0) || (factor < 0.1)) {
        resetWeights = true;
        PRECICE_DEBUG("Significant scaling weight change is detected. The pre-scaling weights will be reset.");
        break;
      }
    }
  }

  if (!resetWeights) {
    return;
  }

  // Reset the weights for non-zero residual sums
  for (size_t k = 0; k < _subVectorSizes.size(); k++) {
    double resSum = _residualSum[k];
    if (not math::equals(resSum, 0.0)) {
      auto offset = _subVectorOffsets[k];
      for (size_t i = 0; i < _subVectorSizes[k]; i++) {
        _weights[i + offset]    = 1 / resSum;
        _invWeights[i + offset] = resSum;
      }
      PRECICE_DEBUG("preconditioner scaling factor[{}] = {}", k, 1 / resSum);
    }
    _previousResidualSum[k] = resSum;
  }
  _requireNewQR = true;
}

} // namespace precice::acceleration::impl
