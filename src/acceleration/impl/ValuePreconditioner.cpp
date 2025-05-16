#include "acceleration/impl/ValuePreconditioner.hpp"
#include <cstddef>
#include <vector>
#include "math/differences.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::acceleration::impl {

ValuePreconditioner::ValuePreconditioner(
    int maxNonConstTimeWindows)
    : Preconditioner(maxNonConstTimeWindows)
{
}

void ValuePreconditioner::_update_(bool                   timeWindowComplete,
                                   const Eigen::VectorXd &oldValues,
                                   const Eigen::VectorXd &res)
{
  if (!timeWindowComplete && !_firstTimeWindow) {
    return;
  }

  std::vector<double> norms(_subVectorSizes.size(), 0.0);

  for (size_t k = 0; k < _subVectorSizes.size(); k++) {
    Eigen::VectorXd part = oldValues.segment(_subVectorOffsets[k], _subVectorSizes[k]);
    norms[k]             = utils::IntraComm::l2norm(part);
  }

  for (size_t k = 0; k < _subVectorSizes.size(); k++) {
    if (norms[k] < math::NUMERICAL_ZERO_DIFFERENCE) {
      PRECICE_WARN("A sub-vector in the residual preconditioner became numerically zero. "
                   "If this occurred in the second iteration and the initial-relaxation factor is equal to 1.0, "
                   "check if the coupling data values of one solver is zero in the first iteration. "
                   "The preconditioner scaling factors were not applied for this iteration.");
    } else {
      for (size_t i = 0; i < _subVectorSizes[k]; i++) {
        PRECICE_ASSERT(norms[k] > 0.0);
        auto offset             = _subVectorOffsets[k];
        _weights[i + offset]    = 1.0 / norms[k];
        _invWeights[i + offset] = norms[k];
      }
    }
  }

  _requireNewQR    = true;
  _firstTimeWindow = false;
}

} // namespace precice::acceleration::impl
