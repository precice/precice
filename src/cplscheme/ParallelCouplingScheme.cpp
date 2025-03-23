#include "ParallelCouplingScheme.hpp"

#include <utility>

#include "cplscheme/BiCouplingScheme.hpp"
#include "logging/LogMacros.hpp"

namespace precice::cplscheme {

ParallelCouplingScheme::ParallelCouplingScheme(
    double             maxTime,
    int                maxTimeWindows,
    double             timeWindowSize,
    const std::string &firstParticipant,
    const std::string &secondParticipant,
    const std::string &localParticipant,
    m2n::PtrM2N        m2n,
    CouplingMode       cplMode,
    int                minIterations,
    int                maxIterations)
    : BiCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, firstParticipant,
                       secondParticipant, localParticipant, std::move(m2n), minIterations, maxIterations, cplMode, constants::TimesteppingMethod::FIXED_TIME_WINDOW_SIZE) {}

ParallelCouplingScheme::ParallelCouplingScheme(
    double             maxTime,
    int                maxTimeWindows,
    double             timeWindowSize,
    const std::string &firstParticipant,
    const std::string &secondParticipant,
    const std::string &localParticipant,
    m2n::PtrM2N        m2n,
    CouplingMode       cplMode)
    : ParallelCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, firstParticipant, secondParticipant, localParticipant, std::move(m2n), cplMode, UNDEFINED_MAX_ITERATIONS, UNDEFINED_MAX_ITERATIONS){};

void ParallelCouplingScheme::exchangeInitialData()
{
  // F: send, receive, S: receive, send
  PRECICE_ASSERT(math::equals(getTime(), getWindowStartTime()), getTime(), getWindowStartTime());
  if (doesFirstStep()) {
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData());
    }
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData());
      notifyDataHasBeenReceived();
    } else {
      initializeWithZeroInitialData(getReceiveData());
    }
  } else { // second participant
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData());
      notifyDataHasBeenReceived();
    } else {
      initializeWithZeroInitialData(getReceiveData());
    }
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData());
    }
  }
}

// Exchanges data before we actually perform the repartitining
// We also set a flag indicating that this data doesn't need to be exchanged
// when advancing the coupling scheme
void ParallelCouplingScheme::exchangeDirectAccessData()
{
  // F: send, receive, S: receive, send
  PRECICE_ASSERT(math::equals(getTime(), getWindowEndTime()), getTime(), getWindowEndTime());

  // get send and receive map
  auto directSend    = filterDataMap(getSendData(), [](const auto &cplData) { return cplData->isDirectAccessWrittenData; });
  auto directReceive = filterDataMap(getReceiveData(), [](const auto &cplData) { return cplData->isDirectAccessWrittenData; });

  if (doesFirstStep()) {
    if (!directSend.empty()) {
      sendData(getM2N(), directSend);
    }
    if (!directReceive.empty()) {
      receiveData(getM2N(), directReceive);
      // notifyDataHasBeenReceived();
    }
    // else {
    //   initializeWithZeroInitialData(getReceiveData());
    // }
  } else { // second participant
    if (!directReceive.empty()) {
      receiveData(getM2N(), directReceive);
      // notifyDataHasBeenReceived();
    }
    // else {
    //   initializeWithZeroInitialData(getReceiveData());
    // }
    if (!directSend.empty()) {
      sendData(getM2N(), directSend);
    }
  }

  _directAccessDataWasSent = true;
}

void ParallelCouplingScheme::exchangeFirstData()
{
  PRECICE_ASSERT(math::equals(getTime(), getWindowEndTime()), getTime(), getWindowEndTime());

  auto relevantSendData    = filterDataMap(getSendData(), [&](const auto &cplData) {
    return !_directAccessDataWasSent || !cplData->isDirectAccessWrittenData;
  });
  auto relevantReceiveData = filterDataMap(getReceiveData(), [&](const auto &cplData) {
    return !_directAccessDataWasSent || !cplData->isDirectAccessWrittenData;
  });

  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), relevantSendData);
  } else { // second participant
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), relevantReceiveData);
    notifyDataHasBeenReceived();
  }
}

void ParallelCouplingScheme::exchangeSecondData()
{
  PRECICE_ASSERT(math::equals(getTime(), getWindowEndTime()), getTime(), getWindowEndTime());
  auto relevantSendData    = filterDataMap(getSendData(), [&](const auto &cplData) {
    return !_directAccessDataWasSent || !cplData->isDirectAccessWrittenData;
  });
  auto relevantReceiveData = filterDataMap(getReceiveData(), [&](const auto &cplData) {
    return !_directAccessDataWasSent || !cplData->isDirectAccessWrittenData;
  });
  if (isExplicitCouplingScheme()) {
    if (doesFirstStep()) { // first participant
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), relevantReceiveData);
      notifyDataHasBeenReceived();
    } else { // second participant
      PRECICE_DEBUG("Sending data...");
      sendData(getM2N(), relevantSendData);
    }
    moveToNextWindow();
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());
    PRECICE_ASSERT(false, "Not (yet) considered");

    if (doesFirstStep()) { // first participant
      PRECICE_DEBUG("Receiving convergence data...");
      receiveConvergence(getM2N());
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      notifyDataHasBeenReceived();
    } else { // second participant
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      doImplicitStep();
      PRECICE_DEBUG("Sending convergence...");
      sendConvergence(getM2N());
      PRECICE_DEBUG("Sending data...");
      sendData(getM2N(), getSendData());
    }

    if (hasConverged()) {
      moveToNextWindow();
    }

    storeIteration();
  }
  _directAccessDataWasSent = false;
}

DataMap &ParallelCouplingScheme::getAccelerationData()
{
  // ParallelCouplingScheme applies acceleration to all CouplingData
  PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration.");
  return _allData;
}

} // namespace precice::cplscheme
