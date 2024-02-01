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
  if (timeWindowComplete || _firstTimeWindow) {

    std::vector<double> norms(_subVectorSizes.size(), 0.0);

    int offset = 0;
    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      Eigen::VectorXd part = Eigen::VectorXd::Zero(_subVectorSizes[k]);
      for (size_t i = 0; i < _subVectorSizes[k]; i++) {
        part(i) = oldValues(i + offset);
      }
      norms[k] = utils::IntraComm::l2norm(part);
      offset += _subVectorSizes[k];
    }

    offset = 0;
    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      if (norms[k] < math::NUMERICAL_ZERO_DIFFERENCE) {
        PRECICE_WARN("A sub-vector in the residual preconditioner became numerically zero. "
                     "If this occurred in the second iteration and the initial-relaxation factor is equal to 1.0, "
                     "check if the coupling data values of one solver is zero in the first iteration. "
                     "The preconditioner scaling factors were not applied for this iteration.");
      } else {
        for (size_t i = 0; i < _subVectorSizes[k]; i++) {
          PRECICE_ASSERT(norms[k] > 0.0);
          _weights[i + offset]    = 1.0 / norms[k];
          _invWeights[i + offset] = norms[k];
        }
      }
      offset += _subVectorSizes[k];
    }

    _requireNewQR    = true;
    _firstTimeWindow = false;
  }
}

} // namespace precice::acceleration::impl
