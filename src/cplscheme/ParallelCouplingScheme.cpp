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
  // @todo bundle receiveConvergence, specific moveToNextWindow and sendData into one function?
  if (doesFirstStep()) { // first participant
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Receiving convergence data...");
      receiveConvergence(getM2N());
    }
    if (hasConverged() || isExplicitCouplingScheme()) {
      moveToNextWindow();
    }
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  }

  if (isImplicitCouplingScheme()) {
    if (not doesFirstStep()) { // second participant
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      doImplicitStep();
    }
  }

  // @todo bundle sendConvergence, specific moveToNextWindow and sendData into one function?
  if (not doesFirstStep()) {
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Sending convergence...");
      sendConvergence(getM2N());
    }
    if (hasConverged() || isExplicitCouplingScheme()) {
      moveToNextWindow();
    }
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  }

  if (isImplicitCouplingScheme()) {
    storeIteration();
  }
}

const DataMap &ParallelCouplingScheme::getAccelerationData()
{
  // ParallelCouplingScheme applies acceleration to all CouplingData
  PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration.");
  return _allData;
}

} // namespace precice::cplscheme
