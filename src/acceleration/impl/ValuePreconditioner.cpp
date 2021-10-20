#include "acceleration/impl/ValuePreconditioner.hpp"
#include <cstddef>
#include <vector>
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace acceleration {
namespace impl {

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
      norms[k] = utils::MasterSlave::l2norm(part);
      offset += _subVectorSizes[k];
      PRECICE_ASSERT(norms[k] > 0.0);
    }

    offset = 0;
    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      for (size_t i = 0; i < _subVectorSizes[k]; i++) {
        _weights[i + offset]    = 1.0 / norms[k];
        _invWeights[i + offset] = norms[k];
      }
      offset += _subVectorSizes[k];
    }

    _requireNewQR    = true;
    _firstTimeWindow = false;
  }
}

} // namespace impl
} // namespace acceleration
} // namespace precice
