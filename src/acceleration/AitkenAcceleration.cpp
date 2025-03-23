#include "acceleration/AitkenAcceleration.hpp"
#include <Eigen/Core>
#include <boost/range/adaptor/map.hpp>
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
      _primaryDataIDs(std::move(dataIDs)),
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
  // Size for each subvector needed for preconditioner
  std::vector<std::size_t> subVectorSizes;
  // Gather sizes
  std::transform(_primaryDataIDs.cbegin(), _primaryDataIDs.cend(), std::back_inserter(subVectorSizes), [&cplData](const auto &d) { return cplData.at(d)->getSize(); });
  // The accumulated sum
  Eigen::Index entries = std::accumulate(subVectorSizes.cbegin(), subVectorSizes.cend(), static_cast<Eigen::Index>(0));

  // Allocate memory
  _oldResiduals = Eigen::VectorXd::Zero(entries);
  _values       = Eigen::VectorXd::Zero(entries);
  _oldValues    = Eigen::VectorXd::Zero(entries);
  if (_primaryDataIDs.size() > 1) {
    _preconditioner->initialize(subVectorSizes);
  }
}

void AitkenAcceleration::performAcceleration(
    DataMap &cplData,
    double   windowStart,
    double   windowEnd)
{
  PRECICE_TRACE();

  // Compute aitken relaxation factor
  PRECICE_ASSERT(utils::contained(*_primaryDataIDs.begin(), cplData));

  concatenateCouplingData(cplData, _primaryDataIDs, _values, _oldValues, windowStart, windowEnd);

  // Compute current residual = values - oldValues
  Eigen::VectorXd residuals = _values - _oldValues;

  // Compute residual deltas (= residuals - oldResiduals) and store it in _oldResiduals
  Eigen::VectorXd residualDeltas = residuals - _oldResiduals;

  // We need to update the preconditioner in every iteration
  if (_primaryDataIDs.size() > 1) {
    _preconditioner->update(false, _values, residuals);
  }

  // Select/compute aitken factor depending on current iteration count
  if (_iterationCounter == 0) {
    // preconditioner not necessary
    _aitkenFactor = math::sign(_aitkenFactor) * std::min(_initialRelaxation, std::abs(_aitkenFactor));
  } else {
    // If we have more than one data set, we scale the data to get a better approximation
    // of the Aitken factor
    if (_primaryDataIDs.size() > 1) {
      _preconditioner->apply(residualDeltas);
      _preconditioner->apply(_oldResiduals);
    }
    // compute fraction of aitken factor with residuals and residual deltas
    double nominator   = utils::IntraComm::dot(_oldResiduals, residualDeltas);
    double denominator = utils::IntraComm::dot(residualDeltas, residualDeltas);
    _aitkenFactor      = -_aitkenFactor * (nominator / denominator);
  }

  PRECICE_DEBUG("AitkenFactor: {}", _aitkenFactor);

  // Perform relaxation with aitken factor
  applyRelaxation(_aitkenFactor, cplData, windowStart);

  // Store residuals for next iteration
  _oldResiduals = std::move(residuals);

  _iterationCounter++;
}

void AitkenAcceleration::iterationsConverged(
    const DataMap &cplData, double windowStart)
{
  _iterationCounter = 0;
  if (_primaryDataIDs.size() > 1) {
    _preconditioner->update(true, _values, _oldResiduals);
  }
  _oldResiduals = Eigen::VectorXd::Constant(_oldResiduals.size(), std::numeric_limits<double>::max());
}

void AitkenAcceleration::concatenateCouplingData(
    const DataMap &cplData, const std::vector<DataID> &dataIDs, Eigen::VectorXd &targetValues, Eigen::VectorXd &targetOldValues, double windowStart, double windowEnd) const
{
  Eigen::Index offset = 0;

  for (auto id : dataIDs) {
    Eigen::Index size = cplData.at(id)->getSize();

    auto valuesSample    = cplData.at(id)->timeStepsStorage().sample(windowEnd);
    auto oldValuesSample = cplData.at(id)->getPreviousValuesAtTime(windowEnd);

    PRECICE_ASSERT(valuesSample.values().size() == size, valuesSample.values().size(), size);
    PRECICE_ASSERT(oldValuesSample.values().size() == size, oldValuesSample.values().size(), size);
    PRECICE_ASSERT(targetValues.size() >= offset + size, "Target vector was not initialized.", targetValues.size(), offset + size);
    PRECICE_ASSERT(targetOldValues.size() >= offset + size, "Target vector was not initialized.");
    for (Eigen::Index i = 0; i < size; i++) {
      targetValues(i + offset)    = valuesSample.values()(i);
      targetOldValues(i + offset) = oldValuesSample.values()(i);
    }
    offset += size;
  }
}
} // namespace precice::acceleration
