#include "acceleration/AitkenAcceleration.hpp"
#include <Eigen/Core>
#include <limits>
#include <map>
#include <math.h>
#include <memory>
#include <ostream>
#include <stddef.h>
#include "cplscheme/CouplingData.hpp"
#include "logging/LogMacros.hpp"
#include "math/math.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace acceleration {

AitkenAcceleration::AitkenAcceleration(double           initialRelaxation,
                                       std::vector<int> dataIDs)
    : _initialRelaxation(initialRelaxation),
      _dataIDs(dataIDs),
      _aitkenFactor(initialRelaxation)
{
  PRECICE_CHECK((_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
                "Initial relaxation factor for Aitken acceleration has to "
                    << "be larger than zero and smaller or equal to one. Current initial relaxation is: " << _initialRelaxation);
}

void AitkenAcceleration::initialize(DataMap &cplData)
{
  checkDataIDs(cplData);
  size_t entries = 0;
  if (_dataIDs.size() == 1) {
    entries = cplData[_dataIDs.at(0)]->values().size();
  } else {
    PRECICE_ASSERT(_dataIDs.size() == 2);
    entries = cplData[_dataIDs.at(0)]->values().size() +
              cplData[_dataIDs.at(1)]->values().size();
  }
  double          initializer = std::numeric_limits<double>::max();
  Eigen::VectorXd toAppend    = Eigen::VectorXd::Constant(entries, initializer);
  utils::append(_residuals, toAppend);

  // Append column for old values if not done by coupling scheme yet
  for (DataMap::value_type &pair : cplData) {
    int cols = pair.second->oldValues.cols();
    if (cols < 1) {
      PRECICE_ASSERT(pair.second->values().size() > 0, pair.first);
      utils::append(pair.second->oldValues,
                    (Eigen::VectorXd) Eigen::VectorXd::Zero(pair.second->values().size()));
    }
  }
}

void AitkenAcceleration::performAcceleration(
    DataMap &cplData)
{
  PRECICE_TRACE();

  // Compute aitken relaxation factor
  PRECICE_ASSERT(utils::contained(*_dataIDs.begin(), cplData));

  Eigen::VectorXd values;
  Eigen::VectorXd oldValues;
  for (int id : _dataIDs) {
    utils::append(values, cplData[id]->values());
    utils::append(oldValues, (Eigen::VectorXd) cplData[id]->oldValues.col(0));
  }

  // Compute current residuals
  Eigen::VectorXd residuals = values;
  residuals -= oldValues;

  // Compute residual deltas and temporarily store it in _residuals
  Eigen::VectorXd residualDeltas = _residuals;
  residualDeltas *= -1.0;
  residualDeltas += residuals;

  // Select/compute aitken factor depending on current iteration count
  if (_iterationCounter == 0) {
    _aitkenFactor = math::sign(_aitkenFactor) * std::min(_initialRelaxation, std::abs(_aitkenFactor));
  } else {
    // compute fraction of aitken factor with residuals and residual deltas
    double nominator   = utils::MasterSlave::dot(_residuals, residualDeltas);
    double denominator = utils::MasterSlave::dot(residualDeltas, residualDeltas);
    _aitkenFactor      = -_aitkenFactor * (nominator / denominator);
  }

  PRECICE_DEBUG("AitkenFactor: " << _aitkenFactor);

  // Perform relaxation with aitken factor
  double omega         = _aitkenFactor;
  double oneMinusOmega = 1.0 - omega;
  for (DataMap::value_type &pair : cplData) {
    auto &      values    = pair.second->values();
    const auto &oldValues = pair.second->oldValues.col(0);
    values *= omega;
    for (int i = 0; i < values.size(); i++) {
      values(i) += oldValues(i) * oneMinusOmega;
    }
  }

  // Store residuals for next iteration
  _residuals = residuals;

  _iterationCounter++;
}

void AitkenAcceleration::iterationsConverged(
    DataMap &cplData)
{
  _iterationCounter = 0;
  _residuals        = Eigen::VectorXd::Constant(_residuals.size(), std::numeric_limits<double>::max());
}

} // namespace acceleration
} // namespace precice
