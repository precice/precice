#include "ParallelCouplingScheme.hpp"

#include <utility>

#include "cplscheme/BiCouplingScheme.hpp"
#include "logging/LogMacros.hpp"

namespace precice::cplscheme {

ParallelCouplingScheme::ParallelCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    int                           validDigits,
    const std::string &           firstParticipant,
    const std::string &           secondParticipant,
    const std::string &           localParticipant,
    m2n::PtrM2N                   m2n,
    constants::TimesteppingMethod dtMethod,
    CouplingMode                  cplMode,
    int                           maxIterations,
    int                           extrapolationOrder)
    : BiCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
                       secondParticipant, localParticipant, std::move(m2n), maxIterations, cplMode, dtMethod, extrapolationOrder) {}

void ParallelCouplingScheme::performReceiveOfFirstAdvance()
{
  // receive nothing by default do constant extrapolation instead
  for (const DataMap::value_type &pair : getReceiveData()) {
    pair.second->moveTimeStepsStorage();
  }
}

void ParallelCouplingScheme::exchangeFirstData()
{
  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  } else { // second participant
    PRECICE_DEBUG("Receiving data...");
    for (const DataMap::value_type &pair : getReceiveData()) {
      pair.second->clearTimeStepsStorage(true);
    }
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  }
}

void ParallelCouplingScheme::exchangeSecondData()
{
  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Receiving data...");
    if (isImplicitCouplingScheme()) {
      receiveConvergence(getM2N());
    }
    for (const DataMap::value_type &pair : getReceiveData()) {
      pair.second->clearTimeStepsStorage(true);
    }
    receiveData(getM2N(), getReceiveData());
    if (hasConverged()) {
      // received converged result of this window, trigger move
      for (const DataMap::value_type &pair : getReceiveData()) {
        pair.second->moveTimeStepsStorage();
      }
    }
    checkDataHasBeenReceived();
  } else { // second participant
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      doImplicitStep();
      sendConvergence(getM2N());
    }
    if (hasConverged()) {
      // received converged result of this window, trigger move
      for (const DataMap::value_type &pair : getReceiveData()) {
        pair.second->moveTimeStepsStorage();
      }
    }
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  }
}

} // namespace precice::cplscheme
