#include <Eigen/Core>
#include <cmath>
#include <cstddef>
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
#include "time/Time.hpp"
#include "time/Waveform.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {

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
      _eps(std::pow(10.0, -1 * validDigits)),
      _extrapolationOrder(extrapolationOrder)
{
  PRECICE_ASSERT(not((maxTime != UNDEFINED_TIME) && (maxTime < 0.0)),
                 "Maximum time has to be larger than zero.");
  PRECICE_ASSERT(not((maxTimeWindows != UNDEFINED_TIME_WINDOWS) && (maxTimeWindows < 0)),
                 "Maximum number of time windows has to be larger than zero.");
  PRECICE_ASSERT(not((timeWindowSize != UNDEFINED_TIME_WINDOW_SIZE) && (timeWindowSize < 0.0)),
                 "Time window size has to be larger than zero.");
  PRECICE_ASSERT((validDigits >= 1) && (validDigits < 17),
                 "Valid digits of time window size has to be between 1 and 16.");
  if (dtMethod == constants::FIXED_TIME_WINDOW_SIZE) {
    PRECICE_ASSERT(hasTimeWindowSize(),
                   "Time window size has to be given when the fixed time window size method is used.");
  }

  PRECICE_ASSERT((maxIterations > 0) || (maxIterations == time::Time::UNDEFINED_EXTRAPOLATION_ORDER),
                 "Maximal iteration limit has to be larger than zero.");

  if (isExplicitCouplingScheme()) {
    PRECICE_ASSERT(maxIterations == UNDEFINED_MAX_ITERATIONS);
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());
    PRECICE_ASSERT(maxIterations >= 1);
  }

  if (isExplicitCouplingScheme()) {
    PRECICE_ASSERT(_extrapolationOrder == time::Time::UNDEFINED_EXTRAPOLATION_ORDER, "Extrapolation is not allowed for explicit coupling");
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());
    PRECICE_CHECK((_extrapolationOrder == 0) || (_extrapolationOrder == 1) || (_extrapolationOrder == 2),
                  "Extrapolation order has to be  0, 1, or 2.");
  }
}

void BaseCouplingScheme::sendData(const m2n::PtrM2N &m2n, const DataMap &sendData)
{
  PRECICE_TRACE();
  std::vector<int> sentDataIDs;
  PRECICE_ASSERT(m2n.get() != nullptr);
  PRECICE_ASSERT(m2n->isConnected());

  for (const DataMap::value_type &pair : sendData) {
    // Data is actually only send if size>0, which is checked in the derived classes implementaiton
    m2n->send(pair.second->values(), pair.second->getMeshID(), pair.second->getDimensions());

    sentDataIDs.push_back(pair.first);
  }
  PRECICE_DEBUG("Number of sent data sets = {}", sentDataIDs.size());
}

void BaseCouplingScheme::receiveData(const m2n::PtrM2N &m2n, const DataMap &receiveData)
{
  PRECICE_TRACE();
  std::vector<int> receivedDataIDs;
  PRECICE_ASSERT(m2n.get());
  PRECICE_ASSERT(m2n->isConnected());
  for (const DataMap::value_type &pair : receiveData) {
    // Data is only received on ranks with size>0, which is checked in the derived class implementation
    m2n->receive(pair.second->values(), pair.second->getMeshID(), pair.second->getDimensions());

    receivedDataIDs.push_back(pair.first);
  }
  PRECICE_DEBUG("Number of received data sets = {}", receivedDataIDs.size());
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
  _time        = startTime;
  _timeWindows = startTimeWindow;

  if (isImplicitCouplingScheme()) {
    if (not doesFirstStep()) {
      PRECICE_CHECK(not _convergenceMeasures.empty(),
                    "At least one convergence measure has to be defined for "
                    "an implicit coupling scheme.");
      // setup convergence measures
      for (ConvergenceMeasureContext &convergenceMeasure : _convergenceMeasures) {
        DataID dataID = convergenceMeasure.couplingData->getDataID();
        assignDataToConvergenceMeasure(&convergenceMeasure, dataID);
      }
      // reserve memory and initialize data with zero
      initializeStorage();
    }
    requireAction(constants::actionWriteIterationCheckpoint());
    initializeTXTWriters();
  }

  initializeImplementation();

  if (_sendsInitializedData) {
    requireAction(constants::actionWriteInitialData());
  }

  // @todo duplicate code also in BaseCouplingScheme::initializeData().
  if (not _sendsInitializedData && not _receivesInitializedData) {
    if (isImplicitCouplingScheme()) {
      if (not doesFirstStep()) {
        storeDataInWaveforms();
        moveToNextWindow();
      }
    }
  }

  _isInitialized = true;
}

void BaseCouplingScheme::initializeData()
{
  // InitializeData uses the template method pattern (https://en.wikipedia.org/wiki/Template_method_pattern).
  PRECICE_ASSERT(_isInitialized);
  PRECICE_ASSERT(not _initializeDataHasBeenCalled);
  _initializeDataHasBeenCalled = true;
  PRECICE_TRACE("initializeData()");

  if (not _sendsInitializedData && not _receivesInitializedData) {
    PRECICE_INFO("initializeData is skipped since no data has to be initialized.");
    return;
  }

  PRECICE_DEBUG("Initializing Data ...");

  _hasDataBeenReceived = false;

  if (isImplicitCouplingScheme()) {
    storeIteration();
  }

  exchangeInitialData();

  // @todo duplicate code also in BaseCouplingScheme::initialize().
  if (isImplicitCouplingScheme()) {
    if (not doesFirstStep()) {
      storeDataInWaveforms();
      moveToNextWindow();
    }
  }
}

void BaseCouplingScheme::advance()
{
  PRECICE_TRACE(_timeWindows, _time);
  checkCompletenessRequiredActions();
  PRECICE_ASSERT(_isInitialized, "Before calling advance() coupling scheme has to be initialized via initialize().");
  _hasDataBeenReceived  = false;
  _isTimeWindowComplete = false;

  PRECICE_ASSERT(_couplingMode != Undefined);

  if (reachedEndOfTimeWindow()) {

    _timeWindows += 1; // increment window counter. If not converged, will be decremented again later.

    bool convergence = exchangeDataAndAccelerate();

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
      PRECICE_ASSERT(_hasDataBeenReceived);
    }
    _computedTimeWindowPart = 0.0; // reset window
  }
}

void BaseCouplingScheme::storeDataInWaveforms()
{
  PRECICE_TRACE(_timeWindows);
  for (DataMap::value_type &pair : _allData) {
    PRECICE_DEBUG("Store data: {}", pair.first);
    _waveforms[pair.first]->store(pair.second->values());
  }
}

void BaseCouplingScheme::moveToNextWindow()
{
  PRECICE_TRACE(_timeWindows);
  for (DataMap::value_type &pair : getAccelerationData()) {
    PRECICE_DEBUG("Store data: {}", pair.first);
    _waveforms[pair.first]->moveToNextWindow();
    pair.second->values() = _waveforms[pair.first]->getInitialGuess();
  }
}

bool BaseCouplingScheme::hasTimeWindowSize() const
{
  return not math::equals(_timeWindowSize, UNDEFINED_TIME_WINDOW_SIZE);
}

double BaseCouplingScheme::getTimeWindowSize() const
{
  PRECICE_ASSERT(not math::equals(_timeWindowSize, UNDEFINED_TIME_WINDOW_SIZE));
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
  _hasDataBeenReceived = true;
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
  if (not math::equals(_timeWindowSize, UNDEFINED_TIME_WINDOW_SIZE)) {
    remainder = getNextTimestepMaxLength();
  }
  PRECICE_DEBUG("return {}", remainder);
  return remainder;
}

double BaseCouplingScheme::getNextTimestepMaxLength() const
{
  if (math::equals(_timeWindowSize, UNDEFINED_TIME_WINDOW_SIZE)) {
    if (math::equals(_maxTime, UNDEFINED_TIME)) {
      return std::numeric_limits<double>::max();
    } else {
      return _maxTime - _time;
    }
  }
  return _timeWindowSize - _computedTimeWindowPart;
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
  if (_timeWindowSize != UNDEFINED_TIME_WINDOW_SIZE) {
    os << ", time-window-size: " << _timeWindowSize;
  }
  if ((_timeWindowSize != UNDEFINED_TIME_WINDOW_SIZE) || (_maxTime != UNDEFINED_TIME)) {
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

void BaseCouplingScheme::initializeStorage()
{
  PRECICE_TRACE();
  // Reserve storage for all data
  for (DataMap::value_type &pair : _allData) {
    _waveforms[pair.first]->initialize(pair.second->values().size());
    pair.second->storeIteration();
  }
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
  PRECICE_ASSERT(_allData.count(dataID) == 1, "Data with given data ID must exist!");
  convMeasure.couplingData = &(*_allData[dataID]);
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
  if (not utils::MasterSlave::isSlave()) {
    _convergenceWriter->writeData("TimeWindow", _timeWindows - 1);
    _convergenceWriter->writeData("Iteration", _iterations);
  }
  for (const auto &convMeasure : _convergenceMeasures) {
    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    PRECICE_ASSERT(convMeasure.measure.get() != nullptr);

    convMeasure.measure->measure(convMeasure.couplingData->previousIteration(), convMeasure.couplingData->values());

    if (not utils::MasterSlave::isSlave() && convMeasure.doesLogging) {
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
  if (not utils::MasterSlave::isSlave()) {

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
  if (not utils::MasterSlave::isSlave()) {

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
  }
}

void BaseCouplingScheme::determineInitialReceive(BaseCouplingScheme::DataMap &receiveData)
{
  if (anyDataRequiresInitialization(receiveData)) {
    _receivesInitializedData = true;
  }
}

void BaseCouplingScheme::addWaveform(int id, const time::PtrWaveform &ptrWaveform)
{
  WaveformMap::value_type waveformPair = std::make_pair(id, ptrWaveform);
  _waveforms.insert(waveformPair);
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

bool BaseCouplingScheme::doImplicitStep()
{
  storeDataInWaveforms();

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

  return convergence;
}

void BaseCouplingScheme::assignDataToConvergenceMeasure(ConvergenceMeasureContext *convergenceMeasure, DataID dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _allData.find(dataID);
  PRECICE_ASSERT(iter != _allData.end(), "Given data ID does not exist in _allData!");
  convergenceMeasure->couplingData = &(*(iter->second));
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

} // namespace cplscheme
} // namespace precice
