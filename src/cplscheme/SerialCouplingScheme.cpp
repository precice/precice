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

void SerialCouplingScheme::performReceiveOfFirstAdvance()
{
  if (doesFirstStep()) {
    // do nothing
  } else { // second participant
    PRECICE_DEBUG("Receiving mesh changes...");
    auto changes = receiveRemoteChanges();
    PRECICE_CHECK(changes.empty(), "Mesh adaptivity in the first step is forbidden.");

    receiveAndSetTimeWindowSize();
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  }
}

CouplingScheme::ChangedMeshes SerialCouplingScheme::firstSynchronization(const CouplingScheme::ChangedMeshes &changes)
{
  PRECICE_DEBUG("First snyc");
  if (!reachedEndOfTimeWindow()) {
    return {};
  };
  // First synchronization point always sends local changes
  PRECICE_DEBUG("Sending mesh changes...");
  sendLocalChanges(changes);
  return {};
}

void SerialCouplingScheme::exchangeFirstData()
{
  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendTimeWindowSize();
    sendData(getM2N(), getSendData());
  } else { // second participant
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Test Convergence and accelerate...");
      doImplicitStep();
      sendConvergence(getM2N());
    }
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  }
}

CouplingScheme::ChangedMeshes SerialCouplingScheme::secondSynchronization()
{
  PRECICE_DEBUG("Second snyc");
  // Second synchronization point always receives remote changes
  if (!reachedEndOfTimeWindow()) {
    return {};
  };
  if (doesFirstStep() || isCouplingOngoing() || (isImplicitCouplingScheme() && not hasConverged())) {
    PRECICE_DEBUG("Receiving mesh changes...");
    return receiveRemoteChanges();
  } else {
    return {};
  }
}

void SerialCouplingScheme::exchangeSecondData()
{
  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Receiving convergence data...");
    if (isImplicitCouplingScheme()) {
      receiveConvergence(getM2N());
    }
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  } else { // second participant
    // the second participant does not want new data in the last iteration of the last time window
    if (isCouplingOngoing() || (isImplicitCouplingScheme() && not hasConverged())) {
      receiveAndSetTimeWindowSize();
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      checkDataHasBeenReceived();
    }
  }
}

const DataMap SerialCouplingScheme::getAccelerationData()
{
  // SerialCouplingSchemes applies acceleration to send data
  return getSendData();
}

} // namespace precice::cplscheme
