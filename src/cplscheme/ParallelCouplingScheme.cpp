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

void ParallelCouplingScheme::exchangeInitialData()
{
  // F: send, receive, S: receive, send
  if (doesFirstStep()) {
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData());
      sendData(getM2N(), getSendGlobalData());
    }
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData());
      receiveData(getM2N(), getReceiveGlobalData());
      checkDataHasBeenReceived();
    } else {
      initializeWithZeroInitialData(getReceiveData());
      initializeWithZeroInitialData(getReceiveGlobalData());
    }
  } else { // second participant
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData());
      receiveData(getM2N(), getReceiveGlobalData());
      checkDataHasBeenReceived();
    } else {
      initializeWithZeroInitialData(getReceiveData());
      initializeWithZeroInitialData(getReceiveGlobalData());
    }
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData());
      sendData(getM2N(), getSendGlobalData());
    }
  }
}

void ParallelCouplingScheme::exchangeFirstData()
{
  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
    sendData(getM2N(), getSendGlobalData());
  } else { // second participant
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    receiveData(getM2N(), getReceiveGlobalData());
    checkDataHasBeenReceived();
  }
}

void ParallelCouplingScheme::exchangeSecondData()
{
  if (isExplicitCouplingScheme()) {
    moveToNextWindow();

    if (doesFirstStep()) { // first participant
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      receiveData(getM2N(), getReceiveGlobalData());
      checkDataHasBeenReceived();
    } else { // second participant
      PRECICE_DEBUG("Sending data...");
      sendData(getM2N(), getSendData());
      sendData(getM2N(), getSendGlobalData());
    }
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());

    if (doesFirstStep()) { // first participant
      if (isImplicitCouplingScheme()) {
        PRECICE_DEBUG("Receiving convergence data...");
        receiveConvergence(getM2N());
      }
    } else { // second participant
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      doImplicitStep();
      if (isImplicitCouplingScheme()) {
        PRECICE_DEBUG("Sending convergence...");
        sendConvergence(getM2N());
      }
    }

    if (hasConverged()) {
      moveToNextWindow();
    }

    if (doesFirstStep()) { // first participant
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      receiveData(getM2N(), getReceiveGlobalData());
      checkDataHasBeenReceived();
    } else { // second participant
      PRECICE_DEBUG("Sending data...");
      sendData(getM2N(), getSendData());
      sendData(getM2N(), getSendGlobalData());
    }

    storeIteration();
  }
}

const DataMap &ParallelCouplingScheme::getAccelerationData()
{
  // ParallelCouplingScheme applies acceleration to all CouplingData
  PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration.");
  return _allMeshData;
}

} // namespace precice::cplscheme
