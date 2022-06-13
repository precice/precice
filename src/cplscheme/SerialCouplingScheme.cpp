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

namespace precice {
namespace cplscheme {

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
      _participantSetsTimeWindowSize = true;  // not allowed to call setTimeWindowSize anymore.
      PRECICE_ASSERT(not hasTimeWindowSize());
    } else {
      _participantReceivesTimeWindowSize = true;
      PRECICE_ASSERT(not _participantSetsTimeWindowSize);
      PRECICE_ASSERT(not solverSetsTimeWindowSize());
    }
  }
}

bool SerialCouplingScheme::solverSetsTimeWindowSize() const
{
  PRECICE_ASSERT(not(hasTimeWindowSize() && _participantSetsTimeWindowSize));
  return _participantSetsTimeWindowSize;
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

void SerialCouplingScheme::initializeImplementation()
{
  // determine whether initial data needs to be communicated
  determineInitialSend(getSendData());
  determineInitialReceive(getReceiveData());

  // If the second participant initializes data, the first receive for the
  // second participant is done in initializeData() instead of initialize().
  if (not doesFirstStep() && not sendsInitializedData() && isCouplingOngoing()) {
    receiveAndSetTimeWindowSize();
    PRECICE_DEBUG("Receiving data");
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  }
}

void SerialCouplingScheme::exchangeInitialData()
{
  if (doesFirstStep()) {
    PRECICE_ASSERT(not sendsInitializedData(), "First participant cannot send data during initialization.");
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData());
      checkInitialDataHasBeenReceived();
    }
  } else { // second participant
    PRECICE_ASSERT(not receivesInitializedData(), "Only first participant can receive data during initialization.");
    if (sendsInitializedData()) {
      // The second participant sends the initialized data to the first participant
      // here, which receives the data on call of initialize().
      sendData(getM2N(), getSendData());
      receiveAndSetTimeWindowSize();
      // This receive replaces the receive in initialize().
      receiveData(getM2N(), getReceiveData());
      checkDataHasBeenReceived();
    }
  }
}

bool SerialCouplingScheme::exchangeDataAndAccelerate()
{
  bool convergence = true;

  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendTimeWindowSize();
    sendData(getM2N(), getSendData());
    if (isImplicitCouplingScheme()) {
      convergence = receiveConvergence(getM2N());
    }
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  } else { // second participant
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Test Convergence and accelerate...");
      convergence = doImplicitStep();
      sendConvergence(getM2N(), convergence);
    }
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
    // the second participant does not want new data in the last iteration of the last time window
    if (isCouplingOngoing() || (isImplicitCouplingScheme() && not convergence)) {
      receiveAndSetTimeWindowSize();
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      checkDataHasBeenReceived();
    }
  }

  return convergence;
}

} // namespace cplscheme
} // namespace precice
