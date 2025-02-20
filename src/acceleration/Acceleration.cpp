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
    for (auto &time : couplingData.timeStepsStorage().getTimes()) {
      if (!math::greater(time, windowStart)) {
        // skip stamples at beginning of this window or earlier since this is either initial data or already converged data from previous windows
        continue;
      }

      auto sample = couplingData.timeStepsStorage().getSampleAtOrAfter(time);
      auto values = sample.values;
      auto old    = couplingData.getPreviousValuesAtTime(time); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
      values      = values * omega + old.values() * (1.0 - omega);

      if (!couplingData.hasGradient()) {
        couplingData.setSampleAtTime(time, time::Sample(couplingData.getDimensions(), values));
      } else {
        auto gradients    = sample.gradients;
        auto oldGradients = couplingData.getPreviousGradientsAtTime(time); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
        gradients         = gradients * omega + oldGradients * (1.0 - omega);
        couplingData.setSampleAtTime(time, time::Sample(couplingData.getDimensions(), values, gradients));
      }
    }
  }
}
} // namespace precice::acceleration
