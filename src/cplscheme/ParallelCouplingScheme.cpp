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
  return; // no action needed.
}

void ParallelCouplingScheme::exchangeFirstData()
{
  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  } else { // second participant
    PRECICE_DEBUG("Receiving data...");
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
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  } else { // second participant
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      doImplicitStep();
      sendConvergence(getM2N());
    }
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  }
  if (hasConverged()) {
    for (const DataMap::value_type &pair : getAllData()) {
      pair.second->moveTimeStepsStorage();
    }
  }
  if (isImplicitCouplingScheme()) {
    storeIteration();
  }
}

} // namespace precice::cplscheme
