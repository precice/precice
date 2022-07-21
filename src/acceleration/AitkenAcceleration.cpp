#include "acceleration/AitkenAcceleration.hpp"
#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <utility>

#include "cplscheme/CouplingData.hpp"
#include "logging/LogMacros.hpp"
#include "math/math.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::acceleration {

AitkenAcceleration::AitkenAcceleration(double                  initialRelaxation,
                                       std::vector<int>        dataIDs,
                                       impl::PtrPreconditioner preconditioner)
    : _initialRelaxation(initialRelaxation),
      _dataIDs(std::move(dataIDs)),
      _aitkenFactor(initialRelaxation),
      _preconditioner(std::move(preconditioner))
{
  PRECICE_CHECK((_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
                "Initial relaxation factor for Aitken acceleration has to "
                "be larger than zero and smaller or equal to one. "
                "Current initial relaxation is: {}",
                _initialRelaxation);
}

void AitkenAcceleration::initialize(const DataMap &cplData)
{
  checkDataIDs(cplData);
  // Accumulate number of entries
  Eigen::Index entries = std::accumulate(_dataIDs.cbegin(), _dataIDs.cend(), static_cast<Eigen::Index>(0), [&cplData](auto &res, const auto &d) { return res + cplData.at(d)->getSize(); });
  // Allocate memory
  _oldResiduals = Eigen::VectorXd::Zero(entries);
  _values       = Eigen::VectorXd::Zero(entries);
  _oldValues    = Eigen::VectorXd::Zero(entries);
}

void AitkenAcceleration::performAcceleration(
    const DataMap &cplData)
{
  PRECICE_TRACE();

  // Compute aitken relaxation factor
  PRECICE_ASSERT(utils::contained(*_dataIDs.begin(), cplData));

  concatenateCouplingData(cplData, _dataIDs, _values, _oldValues);

  // Compute current residuals
  auto residuals = _values - _oldValues;

  // Compute residual deltas and temporarily store it in _oldResiduals
  Eigen::VectorXd residualDeltas = residuals - _oldResiduals;

  // Select/compute aitken factor depending on current iteration count
  if (_iterationCounter == 0) {
    _aitkenFactor = math::sign(_aitkenFactor) * std::min(_initialRelaxation, std::abs(_aitkenFactor));
  } else {
    // compute fraction of aitken factor with residuals and residual deltas
    double nominator   = utils::IntraComm::dot(_oldResiduals, residualDeltas);
    double denominator = utils::IntraComm::dot(residualDeltas, residualDeltas);
    _aitkenFactor      = -_aitkenFactor * (nominator / denominator);
  }

  PRECICE_DEBUG("AitkenFactor: {}", _aitkenFactor);

  // Perform relaxation with aitken factor
  applyRelaxation(_aitkenFactor, cplData);

  // Store residuals for next iteration
  _oldResiduals = residuals;

  _iterationCounter++;
}

void AitkenAcceleration::iterationsConverged(
    const DataMap &cplData)
{
  _iterationCounter = 0;
  _oldResiduals     = Eigen::VectorXd::Constant(_oldResiduals.size(), std::numeric_limits<double>::max());
}

} // namespace precice::acceleration
