#include "SerialCouplingScheme.hpp"
#include "acceleration/Acceleration.hpp"
#include "m2n/M2N.hpp"
#include "math/math.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/MasterSlave.hpp"

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
    int                           maxIterations)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
                         secondParticipant, localParticipant, m2n, maxIterations, dtMethod)
{
  _couplingMode = cplMode;
  // Coupling mode must be either Explicit or Implicit when using SerialCouplingScheme.
  PRECICE_ASSERT(_couplingMode != Undefined);
  if (_couplingMode == Explicit) {
    PRECICE_ASSERT(maxIterations == 1);
  }
}

void SerialCouplingScheme::initialize(
    double startTime,
    int    startTimeWindow)
{
  PRECICE_TRACE(startTime, startTimeWindow);
  PRECICE_ASSERT(not isInitialized());
  PRECICE_ASSERT(math::greaterEquals(startTime, 0.0), startTime);
  PRECICE_ASSERT(startTimeWindow >= 0, startTimeWindow);
  setTime(startTime);
  setTimeWindows(startTimeWindow);

  if (_couplingMode == Implicit) {
    PRECICE_CHECK(not getSendData().empty(), "No send data configured! Use explicit scheme for one-way coupling.");
    if (not doesFirstStep()) {
      if (not _convergenceMeasures.empty()) {
        setupConvergenceMeasures();       // needs _couplingData configured
        setupDataMatrices(getSendData()); // Reserve memory and initialize data with zero
      }
      if (getAcceleration().get() != nullptr) {
        getAcceleration()->initialize(getSendData()); // Reserve memory, initialize
      }
    } else if (getAcceleration().get() != nullptr && getAcceleration()->getDataIDs().size() > 0) {
      int dataID = *(getAcceleration()->getDataIDs().begin());
      PRECICE_CHECK(getSendData(dataID) == nullptr,
                    "In case of serial coupling, acceleration can be defined for "
                        << "data of second participant only!");
    }
    requireAction(constants::actionWriteIterationCheckpoint());
  }

  for (DataMap::value_type &pair : getSendData()) {
    if (pair.second->initialize) {
      PRECICE_CHECK(not doesFirstStep(), "Only second participant can initialize data!");
      PRECICE_DEBUG("Initialized data to be written");
      setHasToSendInitData(true);
      break;
    }
  }

  for (DataMap::value_type &pair : getReceiveData()) {
    if (pair.second->initialize) {
      PRECICE_CHECK(doesFirstStep(), "Only first participant can receive initial data!");
      PRECICE_DEBUG("Initialized data to be received");
      setHasToReceiveInitData(true);
    }
  }

  // If the second participant initializes data, the first receive for the
  // second participant is done in initializeData() instead of initialize().
  if (not doesFirstStep() && not hasToSendInitData() && isCouplingOngoing()) {
    PRECICE_DEBUG("Receiving data");
    receiveAndSetDt();
    receiveData(getM2N());
    setHasDataBeenExchanged(true);
  }

  if (hasToSendInitData()) {
    requireAction(constants::actionWriteInitialData());
  }

  initializeTXTWriters();
  setIsInitialized(true);
}

void SerialCouplingScheme::initializeData()
{
  PRECICE_TRACE();
  PRECICE_CHECK(isInitialized(), "initializeData() can be called after initialize() only!");

  if (not hasToSendInitData() && not hasToReceiveInitData()) {
    PRECICE_INFO("initializeData is skipped since no data has to be initialized");
    return;
  }

  PRECICE_DEBUG("Initializing Data ...");

  PRECICE_CHECK(not(hasToSendInitData() && isActionRequired(constants::actionWriteInitialData())),
                "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);

  if (hasToReceiveInitData() && isCouplingOngoing()) {
    PRECICE_ASSERT(doesFirstStep());
    PRECICE_DEBUG("Receiving data");
    receiveData(getM2N());
    setHasDataBeenExchanged(true);
  }

  if (hasToSendInitData() && isCouplingOngoing()) {
    PRECICE_ASSERT(not doesFirstStep());
    for (DataMap::value_type &pair : getSendData()) {
      if (pair.second->oldValues.cols() == 0)
        break;
      pair.second->oldValues.col(0) = *pair.second->values;
      // For extrapolation, treat the initial value as old time window value
      utils::shiftSetFirst(pair.second->oldValues, *pair.second->values);
    }

    // The second participant sends the initialized data to the first participant
    // here, which receives the data on call of initialize().
    sendData(getM2N());
    receiveAndSetDt();
    // This receive replaces the receive in initialize().
    receiveData(getM2N());
    setHasDataBeenExchanged(true);
  }

  //in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}

void SerialCouplingScheme::advance()
{
  PRECICE_TRACE(getTimeWindows(), getTime());
#ifndef NDEBUG
  for (const DataMap::value_type &pair : getReceiveData()) {
    Eigen::VectorXd &  values = *pair.second->values;
    int                max    = values.size();
    std::ostringstream stream;
    for (int i = 0; (i < max) && (i < 10); i++) {
      stream << values[i] << " ";
    }
    PRECICE_DEBUG("Begin advance, first New Values: " << stream.str());
  }
#endif
  checkCompletenessRequiredActions();

  PRECICE_CHECK(not hasToReceiveInitData() && not hasToSendInitData(),
                "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsTimeWindowComplete(false);

  if (math::equals(getThisTimeWindowRemainder(), 0.0, _eps)) {
    if (_couplingMode == Explicit) {
      explicitAdvance();
    } else if (_couplingMode == Implicit) {
      implicitAdvance();
    }
  } //subcycling completed
}

void SerialCouplingScheme::explicitAdvance()
{
  setIsTimeWindowComplete(true);
  setTimeWindows(getTimeWindows() + 1);
  PRECICE_DEBUG("Sending data...");
  sendDt();
  sendData(getM2N());

  if (isCouplingOngoing() || doesFirstStep()) {
    PRECICE_DEBUG("Receiving data...");
    receiveAndSetDt();
    receiveData(getM2N());
    setHasDataBeenExchanged(true);
  }
  setComputedTimeWindowPart(0.0);
}

void SerialCouplingScheme::implicitAdvance()
{
  bool convergence = true;

  PRECICE_DEBUG("Computed full length of iteration");
  if (doesFirstStep()) {
    sendDt();
    sendData(getM2N());
    getM2N()->receive(convergence);
    if (convergence) {
      timeWindowCompleted();
    }
    receiveData(getM2N());
    setHasDataBeenExchanged(true);
  } else {

    PRECICE_DEBUG("measure convergence.");
    // measure convergence of the coupling iteration,
    convergence = measureConvergence();
    // Stop, when maximal iteration count (given in config) is reached
    if (maxIterationsReached())
      convergence = true;

    // coupling iteration converged for current time window. Advance in time.
    if (convergence) {
      if (getAcceleration().get() != nullptr) {
        _deletedColumnsPPFiltering = getAcceleration()->getDeletedColumns();
        getAcceleration()->iterationsConverged(getSendData());
      }
      newConvergenceMeasurements();
      timeWindowCompleted();

      // no convergence achieved for the coupling iteration within the current time step
    } else if (getAcceleration().get() != nullptr) {
      getAcceleration()->performAcceleration(getSendData());
    }

    // extrapolate new input data for the solver evaluation in time.
    if (convergence && (getExtrapolationOrder() > 0)) {
      extrapolateData(getSendData()); // Also stores data
    } else {                          // Store data for conv. measurement, acceleration, or extrapolation
      for (DataMap::value_type &pair : getSendData()) {
        if (pair.second->oldValues.size() > 0) {
          pair.second->oldValues.col(0) = *pair.second->values;
        }
      }
      for (DataMap::value_type &pair : getReceiveData()) {
        if (pair.second->oldValues.size() > 0) {
          pair.second->oldValues.col(0) = *pair.second->values;
        }
      }
    }

    getM2N()->send(convergence);

    sendData(getM2N());

    // the second participant does not want new data in the last iteration of the last time window
    if (isCouplingOngoing() || not convergence) {
      receiveAndSetDt();
      receiveData(getM2N());
      setHasDataBeenExchanged(true);
    }
  }

  if (not convergence) {
    PRECICE_DEBUG("No convergence achieved");
    requireAction(constants::actionReadIterationCheckpoint());
  } else {
    PRECICE_DEBUG("Convergence achieved");
    advanceTXTWriters();
  }

  updateTimeAndIterations(convergence);
  setComputedTimeWindowPart(0.0);
}

} // namespace cplscheme
} // namespace precice
