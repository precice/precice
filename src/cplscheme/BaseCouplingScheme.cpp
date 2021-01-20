#include "BaseCouplingScheme.hpp"
#include <Eigen/Core>
#include <limits>
#include <math.h>
#include <sstream>
#include <stddef.h>
#include <utility>
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
#include "utils/EigenHelperFunctions.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {

BaseCouplingScheme::BaseCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    int                           validDigits,
    const std::string &           localParticipant,
    int                           maxIterations,
    CouplingMode                  cplMode,
    constants::TimesteppingMethod dtMethod)
    : _couplingMode(cplMode),
      _maxTime(maxTime),
      _maxTimeWindows(maxTimeWindows),
      _timeWindows(1),
      _timeWindowSize(timeWindowSize),
      _maxIterations(maxIterations),
      _iterations(1),
      _totalIterations(1),
      _validDigits(validDigits),
      _localParticipant(localParticipant),
      _eps(std::pow(10.0, -1 * validDigits))
{
  PRECICE_ASSERT(not((maxTime != UNDEFINED_TIME) && (maxTime < 0.0)),
                 "Maximum time has to be larger than zero.");
  PRECICE_ASSERT(not((maxTimeWindows != UNDEFINED_TIME_WINDOWS) && (maxTimeWindows < 0)),
                 "Maximum number of time windows has to be larger than zero.");
  PRECICE_ASSERT(not((timeWindowSize != UNDEFINED_TIME_WINDOW_SIZE) && (timeWindowSize < 0.0)),
                 "Time window size has to be larger than zero.");
  PRECICE_ASSERT((_validDigits >= 1) && (_validDigits < 17),
                 "Valid digits of time window size has to be between 1 and 16.");
  if (dtMethod == constants::FIXED_TIME_WINDOW_SIZE) {
    PRECICE_ASSERT(hasTimeWindowSize(),
                   "Time window size has to be given when the fixed time window size method is used.");
  }

  PRECICE_ASSERT((maxIterations > 0) || (maxIterations == -1),
                 "Maximal iteration limit has to be larger than zero.");

  if (isExplicitCouplingScheme()) {
    PRECICE_ASSERT(maxIterations == -1);
  } else {
    PRECICE_ASSERT(maxIterations >= 1);
  }
}

void BaseCouplingScheme::sendData(m2n::PtrM2N m2n, DataMap sendData)
{
  PRECICE_TRACE();
  std::vector<int> sentDataIDs;
  PRECICE_ASSERT(m2n.get() != nullptr);
  PRECICE_ASSERT(m2n->isConnected());

  for (const DataMap::value_type &pair : sendData) {
    int size = pair.second->values().size();
    if (size > 0) {
      m2n->send(pair.second->values().data(), size, pair.second->mesh->getID(), pair.second->getDimensions());
    }
    sentDataIDs.push_back(pair.first);
  }
  PRECICE_DEBUG("Number of sent data sets = " << sentDataIDs.size());
}

void BaseCouplingScheme::receiveData(m2n::PtrM2N m2n, DataMap receiveData)
{
  PRECICE_TRACE();
  std::vector<int> receivedDataIDs;
  PRECICE_ASSERT(m2n.get());
  PRECICE_ASSERT(m2n->isConnected());
  for (DataMap::value_type &pair : receiveData) {
    int size = pair.second->values().size();
    if (size > 0) {
      m2n->receive(pair.second->values().data(), size, pair.second->mesh->getID(), pair.second->getDimensions());
    }
    receivedDataIDs.push_back(pair.first);
  }
  PRECICE_DEBUG("Number of received data sets = " << receivedDataIDs.size());
}

void BaseCouplingScheme::store(DataMap data)
{
  for (DataMap::value_type &pair : data) {
    if (pair.second->oldValues.size() > 0) {
      pair.second->oldValues.col(0) = pair.second->values();
    }
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
  _time        = startTime;
  _timeWindows = startTimeWindow;

  if (isImplicitCouplingScheme()) {
    if (not doesFirstStep()) {
      PRECICE_CHECK(not _convergenceMeasures.empty(),
                    "At least one convergence measure has to be defined for "
                        << "an implicit coupling scheme.");
      // merge send and receive data for all pp calls
      mergeData();
      // setup convergence measures
      for (ConvergenceMeasureContext &convergenceMeasure : _convergenceMeasures) {
        int dataID = convergenceMeasure.data->getID();
        assignDataToConvergenceMeasure(&convergenceMeasure, dataID);
      }
      // reserve memory and initialize data with zero
      setupDataMatrices(getAccelerationData());
      if (getAcceleration()) {
        getAcceleration()->initialize(getAccelerationData()); // Reserve memory, initialize
      }
    }
    requireAction(constants::actionWriteIterationCheckpoint());
    initializeTXTWriters();
  }

  initializeImplementation();

  if (sendsInitializedData()) {
    requireAction(constants::actionWriteInitialData());
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

  exchangeInitialData();
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

void BaseCouplingScheme::setExtrapolationOrder(
    int order)
{
  PRECICE_CHECK((order == 0) || (order == 1) || (order == 2),
                "Extrapolation order has to be  0, 1, or 2.");
  _extrapolationOrder = order;
}

void BaseCouplingScheme::updateOldValues(DataMap &dataMap)
{
  if (isImplicitCouplingScheme()) {
    for (DataMap::value_type &pair : dataMap) {
      if (pair.second->oldValues.cols() == 0)
        break;
      pair.second->oldValues.col(0) = pair.second->values();
      // For extrapolation, treat the initial value as old time windows value
      utils::shiftSetFirst(pair.second->oldValues, pair.second->values());
    }
  }
}

void BaseCouplingScheme::extrapolateData(DataMap &data)
{
  PRECICE_TRACE(_timeWindows);
  if ((_extrapolationOrder == 1) || getTimeWindows() == 2) { //timesteps is increased before extrapolate is called
    PRECICE_INFO("Performing first order extrapolation");
    for (DataMap::value_type &pair : data) {
      PRECICE_DEBUG("Extrapolate data: " << pair.first);
      PRECICE_ASSERT(pair.second->oldValues.cols() > 1);
      Eigen::VectorXd &values = pair.second->values();
      values *= 2.0;                           // = 2*x^t
      values -= pair.second->oldValues.col(1); // = 2*x^t - x^(t-1)
      utils::shiftSetFirst(pair.second->oldValues, values);
    }
  } else if (_extrapolationOrder == 2) {
    PRECICE_INFO("Performing second order extrapolation");
    for (DataMap::value_type &pair : data) {
      PRECICE_ASSERT(pair.second->oldValues.cols() > 2);
      Eigen::VectorXd &values     = pair.second->values();
      auto             valuesOld1 = pair.second->oldValues.col(1);
      auto             valuesOld2 = pair.second->oldValues.col(2);

      values *= 2.5;              // = 2.5 x^t
      values -= valuesOld1 * 2.0; // = 2.5x^t - 2x^(t-1)
      values += valuesOld2 * 0.5; // = 2.5x^t - 2x^(t-1) + 0.5x^(t-2)
      utils::shiftSetFirst(pair.second->oldValues, values);
    }
  } else {
    PRECICE_ASSERT(false, "Extrapolation order is invalid.");
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
  PRECICE_CHECK(valid, "The timestep "
                       "length given to preCICE in \"advance\" "
                           << timeToAdd << " exceeds the maximum allowed timestep length " << _timeWindowSize - _computedTimeWindowPart + timeToAdd
                           << " in the remaining of this time window. Did you restrict your timestep length, \"dt = min(precice_dt, dt)\" ?"
                           << " For more information, consult the adapter example in the preCICE documentation.");
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
  PRECICE_DEBUG("return " << remainder);
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
    PRECICE_ERROR("The required actions " << stream.str() << " are not fulfilled. Did you forget to call \"markActionFulfilled\"?");
  }
}

void BaseCouplingScheme::setupDataMatrices(DataMap &data)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Data size: " << data.size());
  // Reserve storage for convergence measurement of send and receive data values
  for (ConvergenceMeasureContext &convMeasure : _convergenceMeasures) {
    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    if (convMeasure.couplingData->oldValues.cols() < 1) {
      utils::append(convMeasure.couplingData->oldValues,
                    (Eigen::MatrixXd) Eigen::MatrixXd::Zero(convMeasure.couplingData->values().size(), 1));
    }
  }
  // Reserve storage for extrapolation of data values
  if (_extrapolationOrder > 0) {
    for (DataMap::value_type &pair : data) {
      int cols = pair.second->oldValues.cols();
      PRECICE_DEBUG("Add cols: " << pair.first << ", cols: " << cols);
      PRECICE_ASSERT(cols <= 1, cols);
      utils::append(pair.second->oldValues,
                    (Eigen::MatrixXd) Eigen::MatrixXd::Zero(pair.second->values().size(), _extrapolationOrder + 1 - cols));
    }
  }
  // Storage reservation for acceleration methods happens in Acceleration::initialize
}

void BaseCouplingScheme::setAcceleration(
    acceleration::PtrAcceleration acceleration)
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
    mesh::PtrData               data,
    bool                        suffices,
    bool                        strict,
    impl::PtrConvergenceMeasure measure,
    bool                        doesLogging)
{
  ConvergenceMeasureContext convMeasure;
  convMeasure.data         = std::move(data);
  convMeasure.couplingData = nullptr;
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
  for (size_t i = 0; i < _convergenceMeasures.size(); i++) {
    ConvergenceMeasureContext &convMeasure = _convergenceMeasures[i];

    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    PRECICE_ASSERT(convMeasure.measure.get() != nullptr);
    const auto &oldValues = convMeasure.couplingData->oldValues.col(0);

    convMeasure.measure->measure(oldValues, convMeasure.couplingData->values());

    if (not utils::MasterSlave::isSlave() && convMeasure.doesLogging) {
      _convergenceWriter->writeData(convMeasure.logHeader(), convMeasure.measure->getNormResidual());
    }

    if (not convMeasure.measure->isConvergence()) {
      allConverged = false;
      if (convMeasure.strict) {
        oneStrict = true;
        PRECICE_CHECK(_iterations < _maxIterations,
                      "The strict convergence measure for data \"" + convMeasure.data->getName() +
                          "\" did not converge within the maximum allowed iterations, which terminates the simulation. "
                          "To avoid this forced termination do not mark the convergence measure as strict.")
      }
    } else if (convMeasure.suffices == true) {
      oneSuffices = true;
    }

    PRECICE_INFO(convMeasure.measure->printState());
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
      if (getAcceleration()) {
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

    if (not doesFirstStep() && getAcceleration()) {
      _iterationsWriter->writeData("QNColumns", getAcceleration()->getLSSystemCols());
      _iterationsWriter->writeData("DeletedQNColumns", getAcceleration()->getDeletedColumns());
      _iterationsWriter->writeData("DroppedQNColumns", getAcceleration()->getDroppedColumns());
    }
  }
}

bool BaseCouplingScheme::reachedEndOfTimeWindow()
{
  return math::equals(getThisTimeWindowRemainder(), 0.0, _eps);
}

bool BaseCouplingScheme::maxIterationsReached()
{
  return _iterations == _maxIterations;
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

bool BaseCouplingScheme::accelerate()
{
  PRECICE_DEBUG("measure convergence of the coupling iteration");
  bool convergence = measureConvergence();
  // Stop, when maximal iteration count (given in config) is reached
  if (maxIterationsReached())
    convergence = true;

  // coupling iteration converged for current time window. Advance in time.
  if (convergence) {
    if (getAcceleration()) {
      getAcceleration()->iterationsConverged(getAccelerationData());
    }
    newConvergenceMeasurements();
    // no convergence achieved for the coupling iteration within the current time window
  } else if (getAcceleration()) {
    getAcceleration()->performAcceleration(getAccelerationData());
  }

  // Store data for conv. measurement, acceleration, or extrapolation
  storeData();

  // extrapolate new input data for the solver evaluation in time.
  if (convergence && (_extrapolationOrder > 0)) {
    extrapolateData(getAccelerationData());
  }

  return convergence;
}
} // namespace cplscheme
} // namespace precice
