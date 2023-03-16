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
    auto        storedTimes  = couplingData->getStoredTimesAscending();

    for (auto time : storedTimes) {
      auto data_value = couplingData->getValuesAtTime(time);
      data_value *= omega;

      data_value += oldValues * (1 - omega);
      // Apply relaxation to all timesteps and store it in the current waveform
      couplingData->storeValuesAtTime(time, data_value, true);
    }

    // Ignore the waveform iterations for the gradient for now
    if (couplingData->hasGradient()) {
      auto &      gradients    = couplingData->gradientValues();
      const auto &oldGradients = couplingData->previousIterationGradients();
      gradients *= omega;
      gradients += oldGradients * (1 - omega);
    }
  }
}
} // namespace precice::acceleration
