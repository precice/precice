#include "acceleration/impl/ResidualSumPreconditioner.hpp"
#include <algorithm>
#include <cmath>
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace acceleration {
namespace impl {

ResidualSumPreconditioner::ResidualSumPreconditioner(
    int maxNonConstTimeWindows)
    : Preconditioner(maxNonConstTimeWindows)
{
}

void ResidualSumPreconditioner::initialize(std::vector<size_t> &svs)
{
  PRECICE_TRACE();
  Preconditioner::initialize(svs);

  _residualSum.resize(_subVectorSizes.size(), 0.0);
  _previousScalingWeights.resize(_subVectorSizes.size(), 1.0);
}

void ResidualSumPreconditioner::_update_(bool                   timeWindowComplete,
                                         const Eigen::VectorXd &oldValues,
                                         const Eigen::VectorXd &res)
{
  if (not timeWindowComplete) {
    std::vector<double> norms(_subVectorSizes.size(), 0.0);

    double sum = 0.0;

    int  offset       = 0;
    bool resetWeights = false; // True if pre-scaling weights must be reset
    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      Eigen::VectorXd part = Eigen::VectorXd::Zero(_subVectorSizes[k]);
      for (size_t i = 0; i < _subVectorSizes[k]; i++) {
        part(i) = res(i + offset);
      }
      norms[k] = utils::MasterSlave::dot(part, part);
      sum += norms[k];
      offset += _subVectorSizes[k];
      norms[k] = std::sqrt(norms[k]);
    }
    sum = std::sqrt(sum);
    if (math::equals(sum, 0.0)) {
      PRECICE_WARN("All residual sub-vectors in the residual-sum preconditioner are numerically zero ( sum = {}). "
                   "This indicates that the data values exchanged between two succesive iterations did not change. "
                   "The simulation may be unstable, e.g. produces NAN values. Please check the data values exchanged "
                   "between the solvers is not identical between iterations. The preconditioner scaling factors were "
                   "not updated in this iteration and the scaling factors determined in the previous iteration were used.",
                   sum);
    }

    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      _residualSum[k] += norms[k] / sum;
      if (math::equals(_residualSum[k], 0.0)) {
        PRECICE_WARN("A sub-vector in the residual-sum preconditioner became numerically zero ( sub-vector = {}). "
                     "If this occured in the second iteration and the initial-relaxation factor is equal to 1.0, "
                     "check if the coupling data values of one solver is zero in the first iteration. "
                     "The preconditioner scaling factors were not updated for this iteration and the scaling factors "
                     "determined in the previous iteration were used.",
                     _residualSum[k]);
      }
    }

    offset = 0;
    // normWeights.resize(_subVectorSizes.size());
    // Chech if the new scaling weights are more or less than 1 order of magnitude from the previous weights
    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      double newScalingWeight = (1 / _residualSum[k]);
      if ((newScalingWeight / _previousScalingWeights[k] > 10) || (newScalingWeight / _previousScalingWeights[k] < 0.1)) {
        resetWeights = true;
        _resetSVD    = true;
        PRECICE_INFO("Resetting pre-scaling weights as the value has increased/decreased by more than 1 order of magnitude");
      }
    }

    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      if (not math::equals(_residualSum[k], 0.0)) {
        // Always adjust pre-scaling weights in the first time window
        if (timeWindowPreconditioner < 1 || resetWeights) {
          for (size_t i = 0; i < _subVectorSizes[k]; i++) {
            _weights[i + offset]    = 1 / _residualSum[k];
            _invWeights[i + offset] = _residualSum[k];
          }
          PRECICE_INFO("preconditioner scaling factor[{}] = {}", k, 1 / _residualSum[k]);
          _previousScalingWeights[k] = 1 / _residualSum[k];
          _requireNewQR              = true;
          _areWeightsUpdated         = true;
        }
      }
      //normWeights[k] = 1 / _residualSum[k];
      PRECICE_INFO("Actual Norm of pre-scaling weights in current iteration: {}", _previousScalingWeights[k]);
      offset += _subVectorSizes[k];
    }
    resetWeights = false;

  } else {
    timeWindowPreconditioner++;
    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      _residualSum[k] = 0.0;
    }
  }
}

} // namespace impl
} // namespace acceleration
} // namespace precice
