#include "acceleration/impl/ResidualPreconditioner.hpp"
#include <stddef.h>
#include <vector>
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace acceleration {
namespace impl {

ResidualPreconditioner::ResidualPreconditioner(int maxNonConstTimesteps)
    : Preconditioner(maxNonConstTimesteps)
{
}

void ResidualPreconditioner::_update_(bool                   timestepComplete,
                                      const Eigen::VectorXd &oldValues,
                                      const Eigen::VectorXd &res)
{
  if (not timestepComplete) {
    std::vector<double> norms(_subVectorSizes.size(), 0.0);

    int offset = 0;
    for (size_t k = 0; k < _subVectorSizes.size(); k++) {
      Eigen::VectorXd part = Eigen::VectorXd::Zero(_subVectorSizes[k]);
      for (size_t i = 0; i < _subVectorSizes[k]; i++) {
        part(i) = res(i + offset);
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

    _requireNewQR = true;
  }
}

} // namespace impl
} // namespace acceleration
} // namespace precice
