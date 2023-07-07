#include "SerialCouplingScheme.hpp"
#include <cmath>
#include <memory>
#include <ostream>
#include <utility>

#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/BiCouplingScheme.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/M2N.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme {

SerialCouplingScheme::SerialCouplingScheme(
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
                       secondParticipant, localParticipant, std::move(m2n), maxIterations, cplMode, dtMethod, extrapolationOrder)
{
  if (dtMethod == constants::FIRST_PARTICIPANT_SETS_TIME_WINDOW_SIZE) {
    if (doesFirstStep()) {
      PRECICE_ASSERT(not _participantReceivesTimeWindowSize);
      setTimeWindowSize(UNDEFINED_TIME_WINDOW_SIZE);
      _participantSetsTimeWindowSize = true; // not allowed to call setTimeWindowSize anymore.
      PRECICE_ASSERT(not hasTimeWindowSize());
    } else {
      _participantReceivesTimeWindowSize = true;
      PRECICE_ASSERT(not _participantSetsTimeWindowSize);
    }
  }
}

void SerialCouplingScheme::setTimeWindowSize(double timeWindowSize)
{
  PRECICE_ASSERT(not _participantSetsTimeWindowSize);
  BaseCouplingScheme::setTimeWindowSize(timeWindowSize);
}

void SerialCouplingScheme::sendTimeWindowSize()
{
  PRECICE_TRACE();
  if (_participantSetsTimeWindowSize) {
    PRECICE_DEBUG("sending time window size of {}", getComputedTimeWindowPart());
    getM2N()->send(getComputedTimeWindowPart());
  }
}

double SerialCouplingScheme::getNormalizedWindowTime() const
{
  if (not _participantSetsTimeWindowSize) {
    const double timeWindowStart        = getWindowStartTime();
    const double timeWindowSize         = getTimeWindowSize();
    const double computedTimeWindowPart = getTime() - timeWindowStart;
    // const double computedTimeWindowPart = getComputedTimeWindowPart();  // @todo make public?
    return computedTimeWindowPart / timeWindowSize;
  } else {
    return time::Storage::WINDOW_END; // participant first method does not support subcycling (yet). See https://github.com/precice/precice/issues/1570
  }
}

void SerialCouplingScheme::receiveAndSetTimeWindowSize()
{
  PRECICE_TRACE();
  if (_participantReceivesTimeWindowSize) {
    double dt = UNDEFINED_TIME_WINDOW_SIZE;
    getM2N()->receive(dt);
    PRECICE_DEBUG("Received time window size of {}.", dt);
    PRECICE_ASSERT(not _participantSetsTimeWindowSize);
    PRECICE_ASSERT(not math::equals(dt, UNDEFINED_TIME_WINDOW_SIZE));
    PRECICE_ASSERT(not doesFirstStep(), "Only second participant can receive time window size.");
    setTimeWindowSize(dt);
  }
}

void SerialCouplingScheme::exchangeInitialData()
{
  // F: send, receive, S: receive, send
  if (doesFirstStep()) {
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData());
      receiveData(getM2N(), getReceiveGlobalData());
      checkDataHasBeenReceived();
    } else {
      initializeWithZeroInitialData(getReceiveData());
      initializeWithZeroInitialData(getReceiveGlobalData());
    }
  } else { // second participant
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData());
      sendData(getM2N(), getSendGlobalData());
    }
    // similar to SerialCouplingScheme::exchangeSecondData()
    receiveAndSetTimeWindowSize();
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    receiveData(getM2N(), getReceiveGlobalData());
    checkDataHasBeenReceived();
  }
}

void SerialCouplingScheme::exchangeFirstData()
{
  if (isExplicitCouplingScheme()) {
    if (doesFirstStep()) { // first participant
      PRECICE_DEBUG("Sending data...");
      sendTimeWindowSize();
      sendData(getM2N(), getSendData());
      sendData(getM2N(), getSendGlobalData());
    } else {              // second participant
      moveToNextWindow(); // do moveToNextWindow already here for second participant in SerialCouplingScheme
      PRECICE_DEBUG("Sending data...");
      sendData(getM2N(), getSendData());
      sendData(getM2N(), getSendGlobalData());
    }
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());

    if (doesFirstStep()) { // first participant
      PRECICE_DEBUG("Sending data...");
      sendTimeWindowSize();
      sendData(getM2N(), getSendData());
      sendData(getM2N(), getSendGlobalData());
    } else { // second participant
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      doImplicitStep();
      PRECICE_DEBUG("Sending convergence...");
      sendConvergence(getM2N());
      if (hasConverged()) {
        moveToNextWindow(); // do moveToNextWindow already here for second participant in SerialCouplingScheme
      }
      PRECICE_DEBUG("Sending data...");
      sendData(getM2N(), getSendData());
      sendData(getM2N(), getSendGlobalData());
    }
  }
}

void SerialCouplingScheme::exchangeSecondData()
{
  if (isExplicitCouplingScheme()) {
    if (doesFirstStep()) { // first participant
      moveToNextWindow();  // extrapolation result for receive data of first is directly overwritten in the call of receiveData below
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      receiveData(getM2N(), getReceiveGlobalData());
      checkDataHasBeenReceived();
    }

    if (not doesFirstStep()) { // second participant
      // the second participant does not want new data in the last iteration of the last time window
      if (isCouplingOngoing()) {
        receiveAndSetTimeWindowSize();
        PRECICE_DEBUG("Receiving data...");
        receiveData(getM2N(), getReceiveData());
        receiveData(getM2N(), getReceiveGlobalData());
        checkDataHasBeenReceived();
      }
    }
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());

    if (doesFirstStep()) { // first participant
      PRECICE_DEBUG("Receiving convergence data...");
      receiveConvergence(getM2N());
      if (hasConverged()) {
        moveToNextWindow(); // extrapolation result for receive data of first is directly overwritten in the call of receiveData below
      }
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      receiveData(getM2N(), getReceiveGlobalData());
      checkDataHasBeenReceived();
    }

    storeIteration();

    if (not doesFirstStep()) { // second participant
      // the second participant does not want new data in the last iteration of the last time window
      if (isCouplingOngoing() || not hasConverged()) {
        receiveAndSetTimeWindowSize();
        PRECICE_DEBUG("Receiving data...");
        receiveData(getM2N(), getReceiveData());
        receiveData(getM2N(), getReceiveGlobalData());
        checkDataHasBeenReceived();
      }
    }
  }
}

const DataMap &SerialCouplingScheme::getAccelerationData()
{
  // SerialCouplingSchemes applies acceleration to send data
  return getSendData();
}

} // namespace precice::cplscheme
