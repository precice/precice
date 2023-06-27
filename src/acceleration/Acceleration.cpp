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
      couplingData->sample() = stample.sample;
      auto &      values     = couplingData->values();
      const auto &oldValues  = couplingData->previousIteration(); // @todo need to move previousIteration inside stample? How do deal with changing time grid from one iteration to the next? (Adaptivity?)
      values *= omega;
      values += oldValues * (1 - omega);
      if (couplingData->hasGradient()) {
        auto &      gradients    = couplingData->gradients();
        const auto &oldGradients = couplingData->previousIterationGradients(); // @todo need to move previousIterationGradients inside stample? How do deal with changing time grid from one iteration to the next? (Adaptivity?)
        gradients *= omega;
        gradients += oldGradients * (1 - omega);
      }
      couplingData->setSampleAtTime(stample.timestamp, couplingData->sample());
    }
  }
}
} // namespace precice::acceleration
