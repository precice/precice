#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <sstream>
#include <utility>

#include "BaseCouplingScheme.hpp"
#include "acceleration/Acceleration.hpp"
#include "cplscheme/Constants.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "io/TXTTableWriter.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/IntraComm.hpp"

namespace precice::cplscheme {

BaseCouplingScheme::BaseCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    int                           validDigits,
    std::string                   localParticipant,
    int                           maxIterations,
    CouplingMode                  cplMode,
    constants::TimesteppingMethod dtMethod,
    int                           extrapolationOrder)
    : _couplingMode(cplMode),
      _maxTime(maxTime),
      _maxTimeWindows(maxTimeWindows),
      _timeWindows(1),
      _timeWindowSize(timeWindowSize),
      _maxIterations(maxIterations),
      _iterations(1),
      _totalIterations(1),
      _localParticipant(std::move(localParticipant)),
      _extrapolationOrder(extrapolationOrder),
      _eps(std::pow(10.0, -1 * validDigits))
{
  PRECICE_ASSERT(not((maxTime != UNDEFINED_TIME) && (maxTime < 0.0)),
                 "Maximum time has to be larger than zero.");
  PRECICE_ASSERT(not((maxTimeWindows != UNDEFINED_TIME_WINDOWS) && (maxTimeWindows < 0)),
                 "Maximum number of time windows has to be larger than zero.");
  PRECICE_ASSERT(not(hasTimeWindowSize() && (timeWindowSize < 0.0)),
                 "Time window size has to be larger than zero.");
  PRECICE_ASSERT((validDigits >= 1) && (validDigits < 17),
                 "Valid digits of time window size has to be between 1 and 16.");
  if (dtMethod == constants::FIXED_TIME_WINDOW_SIZE) {
    PRECICE_ASSERT(hasTimeWindowSize(),
                   "Time window size has to be given when the fixed time window size method is used.");
  }

  PRECICE_ASSERT((maxIterations > 0) || (maxIterations == UNDEFINED_MAX_ITERATIONS),
                 "Maximal iteration limit has to be larger than zero.");

  if (isExplicitCouplingScheme()) {
    PRECICE_ASSERT(maxIterations == UNDEFINED_MAX_ITERATIONS);
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());
    PRECICE_ASSERT(maxIterations >= 1);
  }

  if (isExplicitCouplingScheme()) {
    PRECICE_ASSERT(_extrapolationOrder == UNDEFINED_EXTRAPOLATION_ORDER, "Extrapolation is not allowed for explicit coupling");
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());
    PRECICE_CHECK((_extrapolationOrder == 0) || (_extrapolationOrder == 1),
                  "Extrapolation order has to be  0 or 1.");
  }
}

void BaseCouplingScheme::sendNumberOfTimeSteps(const m2n::PtrM2N &m2n, const int numberOfTimeSteps)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Sending number or time steps {}...", numberOfTimeSteps);
  m2n->send(numberOfTimeSteps);
}

void BaseCouplingScheme::sendTimes(const m2n::PtrM2N &m2n, const Eigen::VectorXd times)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Sending times...");
  for (int i = 0; i < times.size(); i++) {
    m2n->send(times(i));
  }
}

void BaseCouplingScheme::sendData(const m2n::PtrM2N &m2n, const DataMap &sendData)
{
  PRECICE_TRACE();
  std::vector<int> sentDataIDs;
  PRECICE_ASSERT(m2n.get() != nullptr);
  PRECICE_ASSERT(m2n->isConnected());

  for (const DataMap::value_type &pair : sendData) {
    auto timesAscending = pair.second->getStoredTimesAscending();
    sendNumberOfTimeSteps(m2n, timesAscending.size());
    sendTimes(m2n, timesAscending);

    PRECICE_ASSERT(math::equals(timesAscending(timesAscending.size() - 1), time::Storage::WINDOW_END), timesAscending(timesAscending.size() - 1)); // assert that last element is time::Storage::WINDOW_END

    auto serializedSamples = pair.second->getSerialized();
    pair.second->clearTimeStepsStorage(false);

    // Data is actually only send if size>0, which is checked in the derived classes implementation
    m2n->send(serializedSamples, pair.second->getMeshID(), pair.second->getDimensions() * timesAscending.size());

    if (pair.second->hasGradient()) {
      m2n->send(pair.second->gradientValues(), pair.second->getMeshID(), pair.second->getDimensions() * pair.second->meshDimensions());
    }

    sentDataIDs.push_back(pair.first);
  }
  PRECICE_DEBUG("Number of sent data sets = {}", sentDataIDs.size());
}

int BaseCouplingScheme::receiveNumberOfTimeSteps(const m2n::PtrM2N &m2n)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Receiving number of time steps...");
  int numberOfTimeSteps;
  m2n->receive(numberOfTimeSteps);
  return numberOfTimeSteps;
}

Eigen::VectorXd BaseCouplingScheme::receiveTimes(const m2n::PtrM2N &m2n, int nTimeSteps)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Receiving times....");
  auto times = Eigen::VectorXd(nTimeSteps);
  for (int i = 0; i < nTimeSteps; i++) {
    m2n->receive(times(i));
  }
  return times;
}

void BaseCouplingScheme::receiveData(const m2n::PtrM2N &m2n, const DataMap &receiveData)
{
  PRECICE_TRACE();
  std::vector<int> receivedDataIDs;
  PRECICE_ASSERT(m2n.get());
  PRECICE_ASSERT(m2n->isConnected());
  for (const DataMap::value_type &pair : receiveData) {
    int nTimeSteps = receiveNumberOfTimeSteps(m2n);

    PRECICE_ASSERT(nTimeSteps > 0);
    auto timesAscending = receiveTimes(m2n, nTimeSteps);

    auto serializedSamples = Eigen::VectorXd(nTimeSteps * pair.second->getSize());
    // Data is only received on ranks with size>0, which is checked in the derived class implementation
    m2n->receive(serializedSamples, pair.second->getMeshID(), pair.second->getDimensions() * nTimeSteps);

    pair.second->storeFromSerialized(timesAscending, serializedSamples);

    if (pair.second->hasGradient()) {
      m2n->receive(pair.second->gradientValues(), pair.second->getMeshID(), pair.second->getDimensions() * pair.second->meshDimensions());
    }

    receivedDataIDs.push_back(pair.first);
  }
  PRECICE_DEBUG("Number of received data sets = {}", receivedDataIDs.size());
  if (receivedDataIDs.size() > 0) {
    _atLeastOneWasReceived = true;
  }
}

void BaseCouplingScheme::initializeZeroReceiveData(const DataMap &receiveData)
{
  for (const DataMap::value_type &pair : receiveData) {
    auto values = pair.second->values();
    pair.second->storeDataAtTime(values, 0.0);
  }
}

void BaseCouplingScheme::setTimeWindowSize(double timeWindowSize)
{
  _timeWindowSize = timeWindowSize;
}

void BaseCouplingScheme::finalize()
{
  PRECICE_TRACE();
  checkCompletenessRequiredActions();
  PRECICE_ASSERT(_isInitialized, "Called finalize() before initialize().");
}

void BaseCouplingScheme::initialize(double startTime, int startTimeWindow)
{
  // Initialize uses the template method pattern (https://en.wikipedia.org/wiki/Template_method_pattern).
  PRECICE_TRACE(startTime, startTimeWindow);
  PRECICE_ASSERT(not isInitialized());
  PRECICE_ASSERT(math::greaterEquals(startTime, 0.0), startTime);
  PRECICE_ASSERT(startTimeWindow >= 0, startTimeWindow);
  _time                = startTime;
  _timeWindows         = startTimeWindow;
  _hasDataBeenReceived = false;

  if (isImplicitCouplingScheme()) {
    if (not doesFirstStep()) {
      PRECICE_CHECK(not _convergenceMeasures.empty(),
                    "At least one convergence measure has to be defined for "
                    "an implicit coupling scheme.");
      // reserve memory and initialize data with zero
      initializeStorages();
    }
    requireAction(constants::actionWriteIterationCheckpoint());
    initializeTXTWriters();
  }

  if (isImplicitCouplingScheme()) {
    storeIteration();
  }

  if (sendsInitializedData()) {
    storeTimeStepSendData(time::Storage::WINDOW_START);
  }

  exchangeInitialData();
  performReceiveOfFirstAdvance();

  _isInitialized = true;
}

bool BaseCouplingScheme::sendsInitializedData() const
{
  return _sendsInitializedData;
}

void BaseCouplingScheme::advance()
{
  PRECICE_TRACE(_timeWindows, _time);
  checkCompletenessRequiredActions();
  PRECICE_ASSERT(_isInitialized, "Before calling advance() coupling scheme has to be initialized via initialize().");
  _hasDataBeenReceived  = false;
  _isTimeWindowComplete = false;

  PRECICE_ASSERT(_couplingMode != Undefined);

  double relativeDt                 = -1;
  bool   usesFirstParticipantMethod = (_timeWindowSize == -1);

  if (isImplicitCouplingScheme()) {
    // store data from writeDataContext in buffer, when advance is called.
    PRECICE_ASSERT(_computedTimeWindowPart > 0);
    if (not usesFirstParticipantMethod) {
      relativeDt = _computedTimeWindowPart / _timeWindowSize;
      PRECICE_ASSERT(math::smallerEquals(relativeDt, time::Storage::WINDOW_END), relativeDt, _computedTimeWindowPart, _timeWindowSize);
      PRECICE_ASSERT(relativeDt > time::Storage::WINDOW_START, relativeDt, _computedTimeWindowPart, _timeWindowSize);
      storeTimeStepSendData(relativeDt);
    } else {
      // We don't support subcycling here, because this is complicated. Therefore, use same strategy like for explicit coupling and just use a single value at end of window.
      // Possible solution: Don't scale times to [0,1], but leave them as they are. Then we would also allow times > 1. We then have two options:
      // 1) scale the times back later when the time window size is known (to still benefit from the simpler handling, if all times are scaled to [0,1]).
      // 2) generally use times in the interval [0, timeWindowSize]. This makes the implementation probably a bit more complicated, but also more consistent.
      if (reachedEndOfTimeWindow()) {                     // only necessary to trigger at end of time window.
        storeTimeStepSendData(time::Storage::WINDOW_END); // only write data at end of window
      }
    }
  } else {
    // work-around for explicit coupling, because it does not support waveform relaxation.
    if (reachedEndOfTimeWindow()) {                     // only necessary to trigger at end of time window.
      storeTimeStepSendData(time::Storage::WINDOW_END); // only write data at end of window
    }
  }

  if (reachedEndOfTimeWindow()) {

    _timeWindows += 1; // increment window counter. If not converged, will be decremented again later.

    bool convergence = exchangeDataAndAccelerate();
    retreiveTimeStepReceiveDataEndOfWindow(); // might be moved into SolverInterfaceImpl.

    if (isImplicitCouplingScheme()) { // check convergence
      if (not convergence) {          // repeat window
        PRECICE_DEBUG("No convergence achieved");
        requireAction(constants::actionReadIterationCheckpoint());
        // The computed time window part equals the time window size, since the
        // time window remainder is zero. Subtract the time window size and do another
        // coupling iteration.
        PRECICE_ASSERT(math::greater(_computedTimeWindowPart, 0.0));
        _time = _time - _computedTimeWindowPart;
        _timeWindows -= 1;
      } else { // write output, prepare for next window
        PRECICE_DEBUG("Convergence achieved");
        advanceTXTWriters();
        PRECICE_INFO("Time window completed");
        _isTimeWindowComplete = true;
        if (isCouplingOngoing()) {
          PRECICE_DEBUG("Setting require create checkpoint");
          requireAction(constants::actionWriteIterationCheckpoint());
        }
      }
      //update iterations
      _totalIterations++;
      if (not convergence) {
        _iterations++;
      } else {
        _iterations = 1;
      }
    } else {
      PRECICE_INFO("Time window completed");
      _isTimeWindowComplete = true;
    }
    if (isCouplingOngoing()) {
      //PRECICE_ASSERT(_hasDataBeenReceived);  // actually incorrect. Data is not necessarily received, if scheme is only sending.
    }
    _computedTimeWindowPart = 0.0; // reset window
  }
}

bool BaseCouplingScheme::hasTimeWindowSize() const
{
  return not math::equals(_timeWindowSize, UNDEFINED_TIME_WINDOW_SIZE);
}

double BaseCouplingScheme::getTimeWindowSize() const
{
  PRECICE_ASSERT(hasTimeWindowSize());
  return _timeWindowSize;
}

void BaseCouplingScheme::addComputedTime(
    double timeToAdd)
{
  PRECICE_TRACE(timeToAdd, _time);
  PRECICE_ASSERT(isCouplingOngoing(), "Invalid call of addComputedTime() after simulation end.");

  // add time interval that has been computed in the solver to get the correct time remainder
  _computedTimeWindowPart += timeToAdd;
  _time += timeToAdd;

  // Check validness
  bool valid = math::greaterEquals(getThisTimeWindowRemainder(), 0.0, _eps);
  PRECICE_CHECK(valid,
                "The timestep length given to preCICE in \"advance\" {} exceeds the maximum allowed timestep length {} "
                "in the remaining of this time window. "
                "Did you restrict your timestep length, \"dt = min(precice_dt, dt)\"? "
                "For more information, consult the adapter example in the preCICE documentation.",
                timeToAdd, _timeWindowSize - _computedTimeWindowPart + timeToAdd);
}

bool BaseCouplingScheme::willDataBeExchanged(
    double lastSolverTimestepLength) const
{
  PRECICE_TRACE(lastSolverTimestepLength);
  double remainder = getThisTimeWindowRemainder() - lastSolverTimestepLength;
  return not math::greater(remainder, 0.0, _eps);
}

bool BaseCouplingScheme::hasDataBeenReceived() const
{
  return _hasDataBeenReceived;
}

void BaseCouplingScheme::checkDataHasBeenReceived()
{
  PRECICE_ASSERT(not _hasDataBeenReceived, "checkDataHasBeenReceived() may only be called once within one coupling iteration. If this assertion is triggered this probably means that your coupling scheme has a bug.");
  if (_atLeastOneWasReceived) {
    _hasDataBeenReceived   = true;
    _atLeastOneWasReceived = false;
  }
}

double BaseCouplingScheme::getTime() const
{
  return _time;
}

int BaseCouplingScheme::getTimeWindows() const
{
  return _timeWindows;
}

double BaseCouplingScheme::getThisTimeWindowRemainder() const
{
  PRECICE_TRACE();
  double remainder = 0.0;
  if (hasTimeWindowSize()) {
    remainder = getNextTimestepMaxLength();
  }
  PRECICE_DEBUG("return {}", remainder);
  return remainder;
}

double BaseCouplingScheme::getNextTimestepMaxLength() const
{
  if (hasTimeWindowSize()) {
    return _timeWindowSize - _computedTimeWindowPart;
  } else {
    if (math::equals(_maxTime, UNDEFINED_TIME)) {
      return std::numeric_limits<double>::max();
    } else {
      return _maxTime - _time;
    }
  }
}

bool BaseCouplingScheme::isCouplingOngoing() const
{
  bool timeLeft      = math::greater(_maxTime, _time, _eps) || math::equals(_maxTime, UNDEFINED_TIME);
  bool timestepsLeft = (_maxTimeWindows >= _timeWindows) || (_maxTimeWindows == UNDEFINED_TIME_WINDOWS);
  return timeLeft && timestepsLeft;
}

bool BaseCouplingScheme::isTimeWindowComplete() const
{
  return _isTimeWindowComplete;
}

bool BaseCouplingScheme::moveWindowBeforeMapping() const
{
  return false; // by default coupling schemes have to move to the next window after performing the mapping
}

void BaseCouplingScheme::retreiveTimeStepReceiveDataEndOfWindow()
{
  if (hasDataBeenReceived()) {
    // needed to avoid problems with round-off-errors.
    auto times       = getReceiveTimes();
    auto largestTime = times.at(times.size() - 1);
    PRECICE_ASSERT(math::equals(largestTime, time::Storage::WINDOW_END), largestTime);
    retreiveTimeStepReceiveData(largestTime); // might be moved into SolverInterfaceImpl.
  }
}

bool BaseCouplingScheme::isActionRequired(
    const std::string &actionName) const
{
  return _actions.count(actionName) > 0;
}

void BaseCouplingScheme::markActionFulfilled(
    const std::string &actionName)
{
  _actions.erase(actionName);
}

void BaseCouplingScheme::requireAction(
    const std::string &actionName)
{
  _actions.insert(actionName);
}

std::string BaseCouplingScheme::printCouplingState() const
{
  std::ostringstream os;
  os << "iteration: " << _iterations; //_iterations;
  if (_maxIterations != -1) {
    os << " of " << _maxIterations;
  }
  os << ", " << printBasicState(_timeWindows, _time) << ", " << printActionsState();
  return os.str();
}

std::string BaseCouplingScheme::printBasicState(
    int    timeWindows,
    double time) const
{
  std::ostringstream os;
  os << "time-window: " << timeWindows;
  if (_maxTimeWindows != UNDEFINED_TIME_WINDOWS) {
    os << " of " << _maxTimeWindows;
  }
  os << ", time: " << time;
  if (_maxTime != UNDEFINED_TIME) {
    os << " of " << _maxTime;
  }
  if (hasTimeWindowSize()) {
    os << ", time-window-size: " << _timeWindowSize;
  }
  if (hasTimeWindowSize() || (_maxTime != UNDEFINED_TIME)) {
    os << ", max-timestep-length: " << getNextTimestepMaxLength();
  }
  os << ", ongoing: ";
  isCouplingOngoing() ? os << "yes" : os << "no";
  os << ", time-window-complete: ";
  _isTimeWindowComplete ? os << "yes" : os << "no";
  return os.str();
}

std::string BaseCouplingScheme::printActionsState() const
{
  std::ostringstream os;
  for (const std::string &actionName : _actions) {
    os << actionName << ' ';
  }
  return os.str();
}

void BaseCouplingScheme::checkCompletenessRequiredActions()
{
  PRECICE_TRACE();
  if (not _actions.empty()) {
    std::ostringstream stream;
    for (const std::string &action : _actions) {
      if (not stream.str().empty()) {
        stream << ", ";
      }
      stream << action;
    }
    PRECICE_ERROR("The required actions {} are not fulfilled. Did you forget to call \"markActionFulfilled\"?", stream.str());
  }
}

void BaseCouplingScheme::initializeStorages()
{
  PRECICE_TRACE();
  // Reserve storage for acceleration
  if (_acceleration) {
    _acceleration->initialize(getAccelerationData());
  }
}

void BaseCouplingScheme::setAcceleration(
    const acceleration::PtrAcceleration &acceleration)
{
  PRECICE_ASSERT(acceleration.get() != nullptr);
  _acceleration = acceleration;
}

void BaseCouplingScheme::newConvergenceMeasurements()
{
  PRECICE_TRACE();
  for (ConvergenceMeasureContext &convMeasure : _convergenceMeasures) {
    PRECICE_ASSERT(convMeasure.measure.get() != nullptr);
    convMeasure.measure->newMeasurementSeries();
  }
}

void BaseCouplingScheme::addConvergenceMeasure(
    int                         dataID,
    bool                        suffices,
    bool                        strict,
    impl::PtrConvergenceMeasure measure,
    bool                        doesLogging)
{
  ConvergenceMeasureContext convMeasure;
  auto                      allData = getAllData();
  PRECICE_ASSERT(allData.count(dataID) == 1, "Data with given data ID must exist!");
  convMeasure.couplingData = allData.at(dataID);
  convMeasure.suffices     = suffices;
  convMeasure.strict       = strict;
  convMeasure.measure      = std::move(measure);
  convMeasure.doesLogging  = doesLogging;
  _convergenceMeasures.push_back(convMeasure);
}

bool BaseCouplingScheme::measureConvergence()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not doesFirstStep());
  bool allConverged = true;
  bool oneSuffices  = false; //at least one convergence measure suffices and did converge
  bool oneStrict    = false; //at least one convergence measure is strict and did not converge
  PRECICE_ASSERT(_convergenceMeasures.size() > 0);
  if (not utils::IntraComm::isSecondary()) {
    _convergenceWriter->writeData("TimeWindow", _timeWindows - 1);
    _convergenceWriter->writeData("Iteration", _iterations);
  }
  for (const auto &convMeasure : _convergenceMeasures) {
    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    PRECICE_ASSERT(convMeasure.measure.get() != nullptr);

    convMeasure.measure->measure(convMeasure.couplingData->previousIteration(), convMeasure.couplingData->values());

    if (not utils::IntraComm::isSecondary() && convMeasure.doesLogging) {
      _convergenceWriter->writeData(convMeasure.logHeader(), convMeasure.measure->getNormResidual());
    }

    if (not convMeasure.measure->isConvergence()) {
      allConverged = false;
      if (convMeasure.strict) {
        oneStrict = true;
        PRECICE_CHECK(_iterations < _maxIterations,
                      "The strict convergence measure for data \"" + convMeasure.couplingData->getDataName() +
                          "\" did not converge within the maximum allowed iterations, which terminates the simulation. "
                          "To avoid this forced termination do not mark the convergence measure as strict.")
      }
    } else if (convMeasure.suffices == true) {
      oneSuffices = true;
    }

    PRECICE_INFO(convMeasure.measure->printState(convMeasure.couplingData->getDataName()));
  }

  if (allConverged) {
    PRECICE_INFO("All converged");
  } else if (oneSuffices && not oneStrict) { //strict overrules suffices
    PRECICE_INFO("Sufficient measures converged");
  }

  return allConverged || (oneSuffices && not oneStrict);
}

void BaseCouplingScheme::initializeTXTWriters()
{
  if (not utils::IntraComm::isSecondary()) {

    _iterationsWriter = std::make_shared<io::TXTTableWriter>("precice-" + _localParticipant + "-iterations.log");
    if (not doesFirstStep()) {
      _convergenceWriter = std::make_shared<io::TXTTableWriter>("precice-" + _localParticipant + "-convergence.log");
    }

    _iterationsWriter->addData("TimeWindow", io::TXTTableWriter::INT);
    _iterationsWriter->addData("TotalIterations", io::TXTTableWriter::INT);
    _iterationsWriter->addData("Iterations", io::TXTTableWriter::INT);
    _iterationsWriter->addData("Convergence", io::TXTTableWriter::INT);

    if (not doesFirstStep()) {
      _convergenceWriter->addData("TimeWindow", io::TXTTableWriter::INT);
      _convergenceWriter->addData("Iteration", io::TXTTableWriter::INT);
    }

    if (not doesFirstStep()) {
      for (ConvergenceMeasureContext &convMeasure : _convergenceMeasures) {

        if (convMeasure.doesLogging) {
          _convergenceWriter->addData(convMeasure.logHeader(), io::TXTTableWriter::DOUBLE);
        }
      }
      if (_acceleration) {
        _iterationsWriter->addData("QNColumns", io::TXTTableWriter::INT);
        _iterationsWriter->addData("DeletedQNColumns", io::TXTTableWriter::INT);
        _iterationsWriter->addData("DroppedQNColumns", io::TXTTableWriter::INT);
      }
    }
  }
}

void BaseCouplingScheme::advanceTXTWriters()
{
  if (not utils::IntraComm::isSecondary()) {

    _iterationsWriter->writeData("TimeWindow", _timeWindows - 1);
    _iterationsWriter->writeData("TotalIterations", _totalIterations);
    _iterationsWriter->writeData("Iterations", _iterations);
    int converged = _iterations < _maxIterations ? 1 : 0;
    _iterationsWriter->writeData("Convergence", converged);

    if (not doesFirstStep() && _acceleration) {
      _iterationsWriter->writeData("QNColumns", _acceleration->getLSSystemCols());
      _iterationsWriter->writeData("DeletedQNColumns", _acceleration->getDeletedColumns());
      _iterationsWriter->writeData("DroppedQNColumns", _acceleration->getDroppedColumns());
    }
  }
}

bool BaseCouplingScheme::reachedEndOfTimeWindow()
{
  return math::equals(getThisTimeWindowRemainder(), 0.0, _eps);
}

void BaseCouplingScheme::determineInitialSend(BaseCouplingScheme::DataMap &sendData)
{
  if (anyDataRequiresInitialization(sendData)) {
    _sendsInitializedData = true;
    requireAction(constants::actionWriteInitialData());
  }
}

void BaseCouplingScheme::determineInitialReceive(BaseCouplingScheme::DataMap &receiveData)
{
  if (anyDataRequiresInitialization(receiveData)) {
    _receivesInitializedData = true;
  }
}

int BaseCouplingScheme::getExtrapolationOrder()
{
  return _extrapolationOrder;
}

void BaseCouplingScheme::retreiveTimeStepForData(double relativeDt, DataID dataId)
{
  auto allData   = getAllData();
  auto data      = allData[dataId];
  data->values() = data->getDataAtTime(relativeDt);
}

bool BaseCouplingScheme::anyDataRequiresInitialization(BaseCouplingScheme::DataMap &dataMap) const
{
  /// @todo implement this function using https://en.cppreference.com/w/cpp/algorithm/all_any_none_of
  for (DataMap::value_type &pair : dataMap) {
    if (pair.second->requiresInitialization) {
      return true;
    }
  }
  return false;
}

void BaseCouplingScheme::storeTimeStepAccelerationDataEndOfWindow()
{
  for (auto &anAccelerationData : getAccelerationData()) {
    auto theData = anAccelerationData.second->values();
    anAccelerationData.second->overrideDataAtEndWindowTime(theData);
  }
}

std::vector<double> BaseCouplingScheme::getAccelerationTimes()
{
  //@todo Should ensure that all times vectors actually hold the same times (since otherwise we would have to get times individually per data), but for BiCouplingScheme this should be fine.
  auto times = std::vector<double>();
  for (auto &data : getAccelerationData()) {
    auto timesVec = data.second->getStoredTimesAscending();
    PRECICE_ASSERT(timesVec.size() > 0, timesVec.size());
    for (int i = 0; i < timesVec.size(); i++) {
      times.push_back(timesVec(i));
    }
    return times;
  }
  PRECICE_ASSERT(false);
}

void BaseCouplingScheme::retreiveTimeStepAccelerationDataEndOfWindow()
{
  // needed to avoid problems with round-off-errors.
  auto times       = getAccelerationTimes();
  auto largestTime = times.at(times.size() - 1);
  PRECICE_ASSERT(math::equals(largestTime, time::Storage::WINDOW_END), largestTime);
  for (auto &anAccelerationData : getAccelerationData()) {
    retreiveTimeStepForData(largestTime, anAccelerationData.second->getDataID());
  }
}

bool BaseCouplingScheme::doImplicitStep()
{
  retreiveTimeStepAccelerationDataEndOfWindow(); // will be needed by acceleration

  PRECICE_DEBUG("measure convergence of the coupling iteration");
  bool convergence = measureConvergence();
  // Stop, when maximal iteration count (given in config) is reached
  if (_iterations == _maxIterations)
    convergence = true;

  // coupling iteration converged for current time window. Advance in time.
  if (convergence) {
    if (_acceleration) {
      _acceleration->iterationsConverged(getAccelerationData());
    }
    newConvergenceMeasurements();
  } else {
    // no convergence achieved for the coupling iteration within the current time window
    if (_acceleration) {
      _acceleration->performAcceleration(getAccelerationData());
    }
  }

  if (convergence) {
    moveToNextWindow();
  }

  // Store data for conv. measurement, acceleration
  storeIteration();

  // Override data with accelerated data
  storeTimeStepAccelerationDataEndOfWindow();

  return convergence;
}

void BaseCouplingScheme::sendConvergence(const m2n::PtrM2N &m2n, bool convergence)
{
  PRECICE_ASSERT(not doesFirstStep(), "For convergence information the sending participant is never the first one.");
  m2n->send(convergence);
}

bool BaseCouplingScheme::receiveConvergence(const m2n::PtrM2N &m2n)
{
  PRECICE_ASSERT(doesFirstStep(), "For convergence information the receiving participant is always the first one.");
  bool convergence;
  m2n->receive(convergence);
  return convergence;
}

} // namespace precice::cplscheme
