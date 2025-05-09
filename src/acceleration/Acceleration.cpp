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

void Acceleration::applyRelaxation(double omega, DataMap &cplData, double windowStart)
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
      if (!math::greater(stample.timestamp, windowStart)) {
        // skip stamples at beginning of this window or earlier since this is either initial data or already converged data from previous windows
        continue;
      }

      auto &values = stample.sample.values;
      auto  old    = couplingData.getPreviousValuesAtTime(stample.timestamp); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
      values       = values * omega + old.values() * (1.0 - omega);

      if (couplingData.hasGradient()) {
        auto &gradients    = stample.sample.gradients;
        auto  oldGradients = couplingData.getPreviousGradientsAtTime(stample.timestamp); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
        gradients          = gradients * omega + oldGradients * (1.0 - omega);
      }
      couplingData.setGlobalSample(stample.sample); // @todo try to remove, currently we have to set CouplingData::_data::_sample to the final stample in the storage. Should not be needed anymore when https://github.com/precice/precice/issues/1645 is resolved.
    }
  }
}
} // namespace precice::acceleration
