#include <acceleration/Acceleration.hpp>
#include "cplscheme/CouplingData.hpp"
#include "utils/Helpers.hpp"

namespace precice::acceleration {

void Acceleration::checkDataIDs(const DataMap &cplData) const
{
#ifndef NDEBUG
  for (int id : getPrimaryDataIDs()) {
    bool valid = utils::contained(id, cplData);
    PRECICE_ASSERT(valid, "Data with ID {} unknown.", id);
  }
#endif
}

void Acceleration::applyRelaxation(double omega, DataMap &cplData)
{
  for (auto &pair : cplData) {
    auto &couplingData = *(pair.second);

    if (couplingData.timeStepsStorage().empty()) {
      continue;
    }
    // @todo: This will apply the scaling to the sample at t=0 and t=1 when
    // calling performAcceleration the first time. Computations at the
    // previousValuesAtTime/oldValues don't change anything and are unneeded
    for (auto &stample : couplingData.timeStepsStorage().stamples()) {
      auto &values    = stample.sample.values;
      auto  oldValues = couplingData.getPreviousValuesAtTime(stample.timestamp); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
      values          = values * omega + oldValues * (1.0 - omega);

      if (couplingData.hasGradient()) {
        auto &gradients    = stample.sample.gradients;
        auto  oldGradients = couplingData.getPreviousGradientsAtTime(stample.timestamp); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
        gradients          = gradients * omega + oldGradients * (1.0 - omega);
      }
    }

    // @todo remove
    // update the "sample"
    couplingData.sample() = couplingData.timeStepsStorage().last().sample;
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
