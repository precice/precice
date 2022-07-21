#include <acceleration/Acceleration.hpp>
#include "cplscheme/CouplingData.hpp"
#include "utils/Helpers.hpp"

namespace precice::acceleration {

void Acceleration::checkDataIDs(const DataMap &cplData) const
{
#ifndef NDEBUG
  for (int id : getDataIDs()) {
    bool valid = utils::contained(id, cplData);
    PRECICE_ASSERT(valid, "Data with ID {} unknown.", id);
  }
#endif
}

void Acceleration::applyRelaxation(double omega, const DataMap &cplData) const
{
  for (const DataMap::value_type &pair : cplData) {
    const auto  couplingData = pair.second;
    auto &      values       = couplingData->values();
    const auto &oldValues    = couplingData->previousIteration();
    values *= omega;
    values += oldValues * (1 - omega);
    if (couplingData->hasGradient()) {
      auto &      gradients    = couplingData->gradientValues();
      const auto &oldGradients = couplingData->previousIterationGradients();
      gradients *= omega;
      gradients += oldGradients * (1 - omega);
    }
  }
}

void Acceleration::concatenateCouplingData(
    const DataMap &cplData, const std::vector<DataID> &dataIDs, Eigen::VectorXd &targetValues, Eigen::VectorXd &targetOldValues) const
{
  Eigen::Index offset = 0;
  for (auto id : dataIDs) {
    Eigen::Index size      = cplData.at(id)->values().size();
    auto &       values    = cplData.at(id)->values();
    const auto & oldValues = cplData.at(id)->previousIteration();
    PRECICE_ASSERT(targetValues.size() >= offset + size, "Target vector was not initialized.", targetValues.size(), offset + size);
    PRECICE_ASSERT(targetOldValues.size() >= offset + size, "Target vector was not initialized.");
    for (Eigen::Index i = 0; i < size; i++) {
      targetValues(i + offset)    = values(i);
      targetOldValues(i + offset) = oldValues(i);
    }
    offset += size;
  }
}
} // namespace precice::acceleration
