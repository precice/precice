#include <Eigen/Core>
#include <algorithm>
#include <boost/range/adaptor/map.hpp>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <sstream>
#include <utility>

#include "BaseCouplingScheme.hpp"
#include "acceleration/Acceleration.hpp"
#include "com/SerializedStamples.hpp"
#include "cplscheme/Constants.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/GlobalCouplingData.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "io/TXTTableWriter.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "mesh/GlobalData.hpp"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

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
  PRECICE_ASSERT(not((maxTime != UNDEFINED_MAX_TIME) && (maxTime < 0.0)),
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
    _extrapolationOrder = 0; // set extrapolation order to zero.
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());
    PRECICE_CHECK((_extrapolationOrder == 0) || (_extrapolationOrder == 1),
                  "Extrapolation order has to be 0 or 1.");
  }
}

bool BaseCouplingScheme::isImplicitCouplingScheme() const
{
  PRECICE_ASSERT(_couplingMode != Undefined);
  return _couplingMode == Implicit;
}

bool BaseCouplingScheme::hasConverged() const
{
  return _hasConverged;
}

void BaseCouplingScheme::sendNumberOfTimeSteps(const m2n::PtrM2N &m2n, const int numberOfTimeSteps)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Sending number or time steps {}...", numberOfTimeSteps);
  m2n->send(numberOfTimeSteps);
}

void BaseCouplingScheme::sendTimes(const m2n::PtrM2N &m2n, const Eigen::VectorXd &times)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Sending times...");
  m2n->send(times);
}

void BaseCouplingScheme::sendData(const m2n::PtrM2N &m2n, const DataMap &sendData)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(m2n.get() != nullptr);
  PRECICE_ASSERT(m2n->isConnected());

  for (const auto &data : sendData | boost::adaptors::map_values) {
    const auto &stamples = data->stamples();
    PRECICE_ASSERT(stamples.size() > 0);

    if (data->exchangeSubsteps()) {
      int nTimeSteps = data->timeStepsStorage().nTimes();
      PRECICE_ASSERT(nTimeSteps > 0);

      if (nTimeSteps > 1) { // otherwise PRECICE_ASSERT below is triggered during initialization
        const Eigen::VectorXd timesAscending = data->timeStepsStorage().getTimes();
        PRECICE_ASSERT(math::equals(timesAscending(0), time::Storage::WINDOW_START), timesAscending(0));                                     // assert that first element is time::Storage::WINDOW_START
        PRECICE_ASSERT(math::equals(timesAscending(nTimeSteps - 1), time::Storage::WINDOW_END), timesAscending(nTimeSteps - 1), nTimeSteps); // assert that last element is time::Storage::WINDOW_END
        sendNumberOfTimeSteps(m2n, nTimeSteps);
        sendTimes(m2n, timesAscending);
      } else { // needs special treatment during initialization
        PRECICE_ASSERT(nTimeSteps == 1);
        nTimeSteps = 2;
        Eigen::VectorXd timesAscending(nTimeSteps);
        timesAscending << time::Storage::WINDOW_START, time::Storage::WINDOW_END;
        sendNumberOfTimeSteps(m2n, nTimeSteps);
        sendTimes(m2n, timesAscending);
      }

      const auto serialized = com::serialize::SerializedStamples::serialize(data);

      // Data is actually only send if size>0, which is checked in the derived classes implementation
      m2n->send(serialized.values(), data->getMeshID(), data->getDimensions() * serialized.nTimeSteps());

      if (data->hasGradient()) {
        m2n->send(serialized.gradients(), data->getMeshID(), data->getDimensions() * data->meshDimensions() * serialized.nTimeSteps());
      }
    } else {
      data->sample() = stamples.back().sample;

      // Data is only received on ranks with size>0, which is checked in the derived class implementation
      m2n->send(data->values(), data->getMeshID(), data->getDimensions());

      if (data->hasGradient()) {
        PRECICE_ASSERT(data->hasGradient());
        m2n->send(data->gradients(), data->getMeshID(), data->getDimensions() * data->meshDimensions());
      }
    }
  }
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
  Eigen::VectorXd times(nTimeSteps);
  m2n->receive(times);
  PRECICE_DEBUG("Received times {}", times);
  return times;
}

void BaseCouplingScheme::receiveData(const m2n::PtrM2N &m2n, const DataMap &receiveData)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(m2n.get());
  PRECICE_ASSERT(m2n->isConnected());
  for (const auto &data : receiveData | boost::adaptors::map_values) {

    if (data->exchangeSubsteps()) {
      const int nTimeSteps = receiveNumberOfTimeSteps(m2n);

      Eigen::VectorXd serializedValues(nTimeSteps * data->getSize());
      PRECICE_ASSERT(nTimeSteps > 0);
      const Eigen::VectorXd timesAscending = receiveTimes(m2n, nTimeSteps);

      auto serialized = com::serialize::SerializedStamples::empty(timesAscending, data);

      // Data is only received on ranks with size>0, which is checked in the derived class implementation
      m2n->receive(serialized.values(), data->getMeshID(), data->getDimensions() * nTimeSteps);

      if (data->hasGradient()) {
        m2n->receive(serialized.gradients(), data->getMeshID(), data->getDimensions() * data->meshDimensions() * nTimeSteps);
      }

      serialized.deserializeInto(timesAscending, data);
    } else {
      // Data is only received on ranks with size>0, which is checked in the derived class implementation
      m2n->receive(data->values(), data->getMeshID(), data->getDimensions());

      if (data->hasGradient()) {
        PRECICE_ASSERT(data->hasGradient());
        m2n->receive(data->gradients(), data->getMeshID(), data->getDimensions() * data->meshDimensions());
      }
      data->setSampleAtTime(time::Storage::WINDOW_END, data->sample());
    }
  }
}

void BaseCouplingScheme::initializeWithZeroInitialData(const DataMap &receiveData)
{
  for (const auto &data : receiveData | boost::adaptors::map_values) {
    PRECICE_DEBUG("Initialize {} as zero.", data->getDataName());
    // just store already initialized zero sample to storage.
    data->setSampleAtTime(time::Storage::WINDOW_START, data->sample());
    data->setSampleAtTime(time::Storage::WINDOW_END, data->sample());
  }
}

PtrCouplingData BaseCouplingScheme::addCouplingData(const mesh::PtrData &data, mesh::PtrMesh mesh, bool requiresInitialization, bool communicateSubsteps)
{
  int             id = data->getID();
  PtrCouplingData ptrCplData;
  if (!utils::contained(id, _allData)) { // data is not used by this coupling scheme yet, create new CouplingData
    ptrCplData = std::make_shared<CouplingData>(data, std::move(mesh), requiresInitialization, communicateSubsteps, _extrapolationOrder);
    _allData.emplace(id, ptrCplData);
  } else { // data is already used by another exchange of this coupling scheme, use existing CouplingData
    ptrCplData = _allData[id];
  }
  return ptrCplData;
}

PtrGlobalCouplingData BaseCouplingScheme::addGlobalCouplingData(const mesh::PtrGlobalData &data, bool requiresInitialization)
{
  int                   id = data->getID();
  PtrGlobalCouplingData ptrGblCplData;
  if (!utils::contained(id, _allGlobalData)) { // data is not used by this coupling scheme yet, create new GlobalCouplingData
    if (isExplicitCouplingScheme()) {
      ptrGblCplData = std::make_shared<GlobalCouplingData>(data, requiresInitialization);
    } else {
      ptrGblCplData = std::make_shared<GlobalCouplingData>(data, requiresInitialization, getExtrapolationOrder());
    }
    _allGlobalData.emplace(id, ptrGblCplData);
    PRECICE_DEBUG("Added global data {} to _allGlobalData.", data->getName());
  } else { // data is already used by another exchange of this coupling scheme, use existing GlobalCouplingData
    ptrGblCplData = _allGlobalData[id];
  }
  return ptrGblCplData;
}

bool BaseCouplingScheme::isExplicitCouplingScheme()
{
  PRECICE_ASSERT(_couplingMode != Undefined);
  return _couplingMode == Explicit;
}

void BaseCouplingScheme::sendGlobalData(const m2n::PtrM2N &m2n, const GlobalDataMap &sendGlobalData)
{
  PRECICE_TRACE();
  std::vector<int> sentGlobalDataIDs; // for debugging
  PRECICE_ASSERT(m2n.get() != nullptr);
  PRECICE_ASSERT(m2n->isConnected());

  for (const auto &data : sendGlobalData | boost::adaptors::map_values) {
    // Data is actually only send if size>0, which is checked in the derived classes implementation
    m2n->send(data->values(), -1, data->getDimensions()); // TODO meshID=-1 is a makeshift thing here. Fix this.
    sentGlobalDataIDs.push_back(data->getDataID());       // Keep track of sent data (for degbugging)
  }
  PRECICE_DEBUG("Number of sent data sets (global) = {}", sentGlobalDataIDs.size());
}

void BaseCouplingScheme::receiveGlobalData(const m2n::PtrM2N &m2n, const GlobalDataMap &receiveGlobalData)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(m2n.get());
  PRECICE_ASSERT(m2n->isConnected());

  for (const auto &data : receiveGlobalData | boost::adaptors::map_values) {
    // Data is only received on ranks with size>0, which is checked in the derived class implementation
    m2n->receive(data->values(), -1, data->getDimensions()); // TODO meshID=-1 is a makeshift thing here. Fix this.
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
  // initialize with zero data here, might eventually be overwritten in exchangeInitialData
  initializeReceiveDataStorage();
  // Initialize uses the template method pattern (https://en.wikipedia.org/wiki/Template_method_pattern).
  PRECICE_TRACE(startTime, startTimeWindow);
  PRECICE_ASSERT(not isInitialized());
  PRECICE_ASSERT(math::greaterEquals(startTime, 0.0), startTime);
  PRECICE_ASSERT(startTimeWindow >= 0, startTimeWindow);
  _timeWindowStartTime = startTime;
  _timeWindows         = startTimeWindow;
  _hasDataBeenReceived = false;

  if (isImplicitCouplingScheme()) {
    storeIteration();
    if (not doesFirstStep()) {
      PRECICE_CHECK(not _convergenceMeasures.empty(),
                    "At least one convergence measure has to be defined for "
                    "an implicit coupling scheme.");
      // reserve memory and initialize data with zero
      if (_acceleration) {
        _acceleration->initialize(getAccelerationData());
      }
    }
    requireAction(CouplingScheme::Action::WriteCheckpoint);
    initializeTXTWriters();
  }

  exchangeInitialData();

  _isInitialized = true;
}

bool BaseCouplingScheme::sendsInitializedData() const
{
  return _sendsInitializedData;
}

CouplingScheme::ChangedMeshes BaseCouplingScheme::firstSynchronization(const CouplingScheme::ChangedMeshes &changes)
{
  PRECICE_ASSERT(changes.empty());
  return changes;
}

void BaseCouplingScheme::firstExchange()
{
  PRECICE_TRACE(_timeWindows, getTime());
  checkCompletenessRequiredActions();
  PRECICE_ASSERT(_isInitialized, "Before calling advance() coupling scheme has to be initialized via initialize().");
  _hasDataBeenReceived  = false;
  _isTimeWindowComplete = false;

  PRECICE_ASSERT(_couplingMode != Undefined);

  if (reachedEndOfTimeWindow()) {

    _timeWindows += 1; // increment window counter. If not converged, will be decremented again later.

    exchangeFirstData();
  }
}

CouplingScheme::ChangedMeshes BaseCouplingScheme::secondSynchronization()
{
  return {};
}

void BaseCouplingScheme::secondExchange()
{
  PRECICE_TRACE(_timeWindows, getTime());
  checkCompletenessRequiredActions();
  PRECICE_ASSERT(_isInitialized, "Before calling advance() coupling scheme has to be initialized via initialize().");
  PRECICE_ASSERT(_couplingMode != Undefined);

  // from first phase
  PRECICE_ASSERT(!_isTimeWindowComplete);

  if (reachedEndOfTimeWindow()) {

    exchangeSecondData();

    if (isImplicitCouplingScheme()) { // check convergence
      if (not hasConverged()) {       // repeat window
        PRECICE_DEBUG("No convergence achieved");
        requireAction(CouplingScheme::Action::ReadCheckpoint);
        // The computed time window part equals the time window size, since the
        // time window remainder is zero. Subtract the time window size and do another
        // coupling iteration.
        PRECICE_ASSERT(math::greater(_computedTimeWindowPart, 0.0));
        _timeWindows -= 1;
        _computedTimeWindowPart = 0.0; // reset window
      } else {                         // write output, prepare for next window
        PRECICE_DEBUG("Convergence achieved");
        advanceTXTWriters();
        PRECICE_INFO("Time window completed");
        _isTimeWindowComplete = true;
        _timeWindowStartTime += _computedTimeWindowPart;
        _computedTimeWindowPart = 0.0; // reset window
        if (isCouplingOngoing()) {
          PRECICE_DEBUG("Setting require create checkpoint");
          requireAction(CouplingScheme::Action::WriteCheckpoint);
        }
      }
      //update iterations
      _totalIterations++;
      if (not hasConverged()) {
        _iterations++;
      } else {
        _iterations = 1;
      }
    } else {
      PRECICE_INFO("Time window completed");
      _isTimeWindowComplete = true;
      _timeWindowStartTime += _computedTimeWindowPart;
      _computedTimeWindowPart = 0.0; // reset window
    }
    if (isCouplingOngoing()) {
      PRECICE_ASSERT(_hasDataBeenReceived);
    }
  }
}

void BaseCouplingScheme::moveToNextWindow()
{
  PRECICE_TRACE(_timeWindows);
  for (auto &data : _allData | boost::adaptors::map_values) {
    data->moveToNextWindow();
  }
  // TODO: Add global data here?
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

double BaseCouplingScheme::getWindowStartTime() const
{
  return getTime() - _computedTimeWindowPart;
}

double BaseCouplingScheme::getNormalizedWindowTime() const
{
  return getComputedTimeWindowPart() / getTimeWindowSize();
}

bool BaseCouplingScheme::isInitialized() const
{
  return _isInitialized;
}

bool BaseCouplingScheme::addComputedTime(
    double timeToAdd)
{
  PRECICE_TRACE(timeToAdd, getTime());
  PRECICE_ASSERT(isCouplingOngoing(), "Invalid call of addComputedTime() after simulation end.");

  // add time interval that has been computed in the solver to get the correct time remainder
  _computedTimeWindowPart += timeToAdd;

  // Check validness
  bool valid = math::greaterEquals(getNextTimeStepMaxSize(), 0.0, _eps);
  PRECICE_CHECK(valid,
                "The time step size given to preCICE in \"advance\" {} exceeds the maximum allowed time step size {} "
                "in the remaining of this time window. "
                "Did you restrict your time step size, \"dt = min(preciceDt, solverDt)\"? "
                "For more information, consult the adapter example in the preCICE documentation.",
                timeToAdd, _timeWindowSize - _computedTimeWindowPart + timeToAdd);

  const bool isAtWindowEnd = math::equals(getNormalizedWindowTime(), time::Storage::WINDOW_END, _eps);
  return isAtWindowEnd;
}

bool BaseCouplingScheme::willDataBeExchanged(
    double lastSolverTimeStepSize) const
{
  PRECICE_TRACE(lastSolverTimeStepSize);
  double remainder = getNextTimeStepMaxSize() - lastSolverTimeStepSize;
  return not math::greater(remainder, 0.0, _eps);
}

bool BaseCouplingScheme::hasDataBeenReceived() const
{
  return _hasDataBeenReceived;
}

double BaseCouplingScheme::getComputedTimeWindowPart() const
{
  return _computedTimeWindowPart;
}

void BaseCouplingScheme::setDoesFirstStep(bool doesFirstStep)
{
  _doesFirstStep = doesFirstStep;
}

void BaseCouplingScheme::checkDataHasBeenReceived()
{
  PRECICE_ASSERT(not _hasDataBeenReceived, "checkDataHasBeenReceived() may only be called once within one coupling iteration. If this assertion is triggered this probably means that your coupling scheme has a bug.");
  _hasDataBeenReceived = true;
}

bool BaseCouplingScheme::receivesInitializedData() const
{
  return _receivesInitializedData;
}

void BaseCouplingScheme::setTimeWindows(int timeWindows)
{
  _timeWindows = timeWindows;
}

double BaseCouplingScheme::getTime() const
{
  return _timeWindowStartTime + _computedTimeWindowPart;
}

int BaseCouplingScheme::getTimeWindows() const
{
  return _timeWindows;
}

double BaseCouplingScheme::getNextTimeStepMaxSize() const
{
  if (hasTimeWindowSize()) {
    return _timeWindowSize - _computedTimeWindowPart;
  } else {
    if (math::equals(_maxTime, UNDEFINED_MAX_TIME)) {
      return std::numeric_limits<double>::max();
    } else {
      return _maxTime - getTime();
    }
  }
}

bool BaseCouplingScheme::isCouplingOngoing() const
{
  bool timeLeft      = math::greater(_maxTime, getTime(), _eps) || math::equals(_maxTime, UNDEFINED_MAX_TIME);
  bool timestepsLeft = (_maxTimeWindows >= _timeWindows) || (_maxTimeWindows == UNDEFINED_TIME_WINDOWS);
  return timeLeft && timestepsLeft;
}

bool BaseCouplingScheme::isTimeWindowComplete() const
{
  return _isTimeWindowComplete;
}

bool BaseCouplingScheme::isActionRequired(
    Action action) const
{
  return _requiredActions.count(action) == 1;
}

bool BaseCouplingScheme::isActionFulfilled(
    Action action) const
{
  return _fulfilledActions.count(action) == 1;
}

void BaseCouplingScheme::markActionFulfilled(
    Action action)
{
  PRECICE_ASSERT(isActionRequired(action));
  _fulfilledActions.insert(action);
}

void BaseCouplingScheme::requireAction(
    Action action)
{
  _requiredActions.insert(action);
}

std::string BaseCouplingScheme::printCouplingState() const
{
  std::ostringstream os;
  os << "iteration: " << _iterations; //_iterations;
  if (_maxIterations != -1) {
    os << " of " << _maxIterations;
  }
  os << ", " << printBasicState(_timeWindows, getTime()) << ", " << printActionsState();
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
  if (_maxTime != UNDEFINED_MAX_TIME) {
    os << " of " << _maxTime;
  }
  if (hasTimeWindowSize()) {
    os << ", time-window-size: " << _timeWindowSize;
  }
  if (hasTimeWindowSize() || (_maxTime != UNDEFINED_MAX_TIME)) {
    os << ", max-time-step-size: " << getNextTimeStepMaxSize();
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
  for (auto action : _requiredActions) {
    os << toString(action) << ' ';
  }
  return os.str();
}

void BaseCouplingScheme::checkCompletenessRequiredActions()
{
  PRECICE_TRACE();
  std::vector<Action> missing;
  std::set_difference(_requiredActions.begin(), _requiredActions.end(),
                      _fulfilledActions.begin(), _fulfilledActions.end(),
                      std::back_inserter(missing));
  if (not missing.empty()) {
    std::ostringstream stream;
    for (auto action : missing) {
      if (not stream.str().empty()) {
        stream << ", ";
      }
      stream << toString(action);
    }
    PRECICE_ERROR("The required actions {} are not fulfilled. "
                  "Did you forget to call \"requiresReadingCheckpoint()\" or \"requiresWritingCheckpoint()\"?",
                  stream.str());
  }
  _requiredActions.clear();
  _fulfilledActions.clear();
}

void BaseCouplingScheme::setAcceleration(
    const acceleration::PtrAcceleration &acceleration)
{
  PRECICE_ASSERT(acceleration.get() != nullptr);
  _acceleration = acceleration;
}

bool BaseCouplingScheme::doesFirstStep() const
{
  return _doesFirstStep;
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
  convMeasure.couplingData = _allData.at(dataID);
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
  bool oneSuffices  = false; // at least one convergence measure suffices and did converge
  bool oneStrict    = false; // at least one convergence measure is strict and did not converge
  PRECICE_ASSERT(_convergenceMeasures.size() > 0);
  if (not utils::IntraComm::isSecondary()) {
    _convergenceWriter->writeData("TimeWindow", _timeWindows - 1);
    _convergenceWriter->writeData("Iteration", _iterations);
  }
  for (const auto &convMeasure : _convergenceMeasures) {
    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    PRECICE_ASSERT(convMeasure.measure.get() != nullptr);
    PRECICE_ASSERT(convMeasure.couplingData->previousIteration().size() == convMeasure.couplingData->values().size(), convMeasure.couplingData->previousIteration().size(), convMeasure.couplingData->values().size(), convMeasure.couplingData->getDataName());
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
  } else if (oneSuffices && not oneStrict) { // strict overrules suffices
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
  return math::equals(getNextTimeStepMaxSize(), 0.0, _eps) || not hasTimeWindowSize();
}

void BaseCouplingScheme::storeIteration()
{
  PRECICE_ASSERT(isImplicitCouplingScheme());
  for (const auto &data : _allData | boost::adaptors::map_values) {
    data->storeIteration();
  }
  for (const auto &data : _allGlobalData | boost::adaptors::map_values) {
    data->storeIteration();
  }
}

void BaseCouplingScheme::determineInitialSend(DataMap &sendData)
{
  if (anyDataRequiresInitialization(sendData)) {
    _sendsInitializedData = true;
    requireAction(CouplingScheme::Action::InitializeData);
  }
}

void BaseCouplingScheme::determineInitialSend(GlobalDataMap &sendGlobalData)
{
  if (anyDataRequiresInitialization(sendGlobalData)) {
    _sendsInitializedData = true;
    requireAction(CouplingScheme::Action::InitializeData);
  }
  //TODO test this
}

void BaseCouplingScheme::determineInitialReceive(DataMap &receiveData)
{
  if (anyDataRequiresInitialization(receiveData)) {
    _receivesInitializedData = true;
  }
}

bool BaseCouplingScheme::anyDataRequiresInitialization(DataMap &dataMap) const
{
  /// @todo implement this function using https://en.cppreference.com/w/cpp/algorithm/all_any_none_of
  for (const auto &data : dataMap | boost::adaptors::map_values) {
    if (data->requiresInitialization) {
      return true;
    }
  }
  return false;
}

bool BaseCouplingScheme::anyDataRequiresInitialization(GlobalDataMap &globalDataMap) const
{
  /// @todo implement this function using https://en.cppreference.com/w/cpp/algorithm/all_any_none_of
  for (const auto &data : globalDataMap | boost::adaptors::map_values) {
    if (data->requiresInitialization) {
      PRECICE_ERROR("Initialization for global data exchange is not tested yet.");
      return true;
    }
  }
  return false;
  //TODO test this
}

void BaseCouplingScheme::doImplicitStep()
{
  PRECICE_DEBUG("measure convergence of the coupling iteration");
  _hasConverged = measureConvergence();
  // Stop, when maximal iteration count (given in config) is reached
  if (_iterations == _maxIterations)
    _hasConverged = true;

  // coupling iteration converged for current time window. Advance in time.
  if (_hasConverged) {
    if (_acceleration) {
      _acceleration->iterationsConverged(getAccelerationData());
    }
    newConvergenceMeasurements();
  } else {
    // no convergence achieved for the coupling iteration within the current time window
    if (_acceleration) {
      // Acceleration works on CouplingData::values(), so we retrieve the data from the storage, perform the acceleration and then put the data back into the storage. See also https://github.com/precice/precice/issues/1645.
      // @todo For acceleration schemes as described in "Rüth, B, Uekermann, B, Mehl, M, Birken, P, Monge, A, Bungartz, H-J. Quasi-Newton waveform iteration for partitioned surface-coupled multiphysics applications. https://doi.org/10.1002/nme.6443" we need a more elaborate implementation.

      // Load from storage into buffer
      for (auto &data : getAccelerationData() | boost::adaptors::map_values) {
        const auto &stamples = data->stamples();
        PRECICE_ASSERT(stamples.size() > 0);
        data->sample() = stamples.back().sample;
      }

      _acceleration->performAcceleration(getAccelerationData());

      // Store from buffer
      // @todo Currently only data at time::Storage::WINDOW_END is accelerated. Remaining data in storage stays as it is.
      for (auto &data : getAccelerationData() | boost::adaptors::map_values) {
        data->setSampleAtTime(time::Storage::WINDOW_END, data->sample());
      }
    }
  }
}

void BaseCouplingScheme::sendConvergence(const m2n::PtrM2N &m2n)
{
  PRECICE_ASSERT(isImplicitCouplingScheme());
  PRECICE_ASSERT(not doesFirstStep(), "For convergence information the sending participant is never the first one.");
  m2n->send(_hasConverged);
}

void BaseCouplingScheme::receiveConvergence(const m2n::PtrM2N &m2n)
{
  PRECICE_ASSERT(isImplicitCouplingScheme());
  PRECICE_ASSERT(doesFirstStep(), "For convergence information the receiving participant is always the first one.");
  m2n->receive(_hasConverged);
}

} // namespace precice::cplscheme
