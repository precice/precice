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
    int                           maxIterations)
    : BiCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
                       secondParticipant, localParticipant, std::move(m2n), maxIterations, cplMode, dtMethod) {}

void ParallelCouplingScheme::exchangeInitialData()
{
  bool initialCommunication = true;

  // F: send, receive, S: receive, send
  if (doesFirstStep()) {
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData(), initialCommunication);
    }
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData(), initialCommunication);
      notifyDataHasBeenReceived();
    } else {
      initializeWithZeroInitialData(getReceiveData());
    }
  } else { // second participant
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData(), initialCommunication);
      notifyDataHasBeenReceived();
    } else {
      initializeWithZeroInitialData(getReceiveData());
    }
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData(), initialCommunication);
    }
  }
}

void ParallelCouplingScheme::exchangeFirstData()
{
  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  } else { // second participant
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    notifyDataHasBeenReceived();
  }
}

void ParallelCouplingScheme::exchangeSecondData()
{
  if (isExplicitCouplingScheme()) {
    moveToNextWindow();

    if (doesFirstStep()) { // first participant
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      notifyDataHasBeenReceived();
    } else { // second participant
      PRECICE_DEBUG("Sending data...");
      sendData(getM2N(), getSendData());
    }
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());

    if (doesFirstStep()) { // first participant
      PRECICE_DEBUG("Receiving convergence data...");
      receiveConvergence(getM2N());
      if (hasConverged()) {
        moveToNextWindow();
      }
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      notifyDataHasBeenReceived();
    } else { // second participant
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      doImplicitStep();
      PRECICE_DEBUG("Sending convergence...");
      sendConvergence(getM2N());
      if (hasConverged()) {
        moveToNextWindow();
      }
      PRECICE_DEBUG("Sending data...");
      sendData(getM2N(), getSendData());
    }

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
