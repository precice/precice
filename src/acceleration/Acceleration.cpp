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
    const auto couplingData = pair.second;

    for (auto &stample : couplingData->stamples()) {
      auto values            = stample.sample.values;
      auto oldValues         = couplingData->getPreviousValuesAtTime(stample.timestamp); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
      couplingData->values() = values * omega;
      couplingData->values() += oldValues * (1 - omega);

      if (couplingData->hasGradient()) {
        auto       gradients      = stample.sample.gradients;
        const auto oldGradients   = couplingData->getPreviousGradientsAtTime(stample.timestamp); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
        couplingData->gradients() = gradients * omega;
        couplingData->gradients() += oldGradients * (1 - omega);
      }

      couplingData->setSampleAtTime(stample.timestamp, couplingData->sample());
    }
  }
}
} // namespace precice::acceleration
