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

#include "acceleration/Acceleration.hpp"
#include "acceleration/BaseQNAcceleration.hpp"
#include "com/SerializedStamples.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
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
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme {

BaseCouplingScheme::BaseCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    std::string                   localParticipant,
    int                           minIterations,
    int                           maxIterations,
    CouplingMode                  cplMode,
    constants::TimesteppingMethod dtMethod)
    : _couplingMode(cplMode),
      _maxTime(maxTime),
      _maxTimeWindows(maxTimeWindows),
      _timeWindows(1),
      _timeWindowSize(timeWindowSize),
      _nextTimeWindowSize(timeWindowSize),
      _minIterations(minIterations),
      _maxIterations(maxIterations),
      _iterations(1),
      _totalIterations(1),
      _localParticipant(std::move(localParticipant))
{
  PRECICE_ASSERT(not((maxTime != UNDEFINED_MAX_TIME) && (maxTime < 0.0)),
                 "Maximum time has to be larger than zero.");
  PRECICE_ASSERT(not((maxTimeWindows != UNDEFINED_TIME_WINDOWS) && (maxTimeWindows < 0)),
                 "Maximum number of time windows has to be larger than zero.");
  PRECICE_ASSERT(not(hasTimeWindowSize() && (timeWindowSize < 0.0)),
                 "Time window size has to be larger than zero.");
  if (dtMethod == constants::FIXED_TIME_WINDOW_SIZE) {
    PRECICE_ASSERT(hasTimeWindowSize(),
                   "Time window size has to be given when the fixed time window size method is used.");
  }

  if (isExplicitCouplingScheme()) {
    PRECICE_ASSERT(minIterations == UNDEFINED_MIN_ITERATIONS);
    PRECICE_ASSERT(maxIterations == UNDEFINED_MAX_ITERATIONS);
  } else {
    PRECICE_ASSERT(isImplicitCouplingScheme());
    PRECICE_ASSERT(minIterations != UNDEFINED_MIN_ITERATIONS);
    PRECICE_ASSERT(maxIterations != UNDEFINED_MAX_ITERATIONS);

    PRECICE_ASSERT(minIterations > 0,
                   minIterations,
                   "Minimal iteration limit has to be larger than zero.");
    PRECICE_ASSERT((maxIterations == INFINITE_MAX_ITERATIONS) || (maxIterations > 0),
                   maxIterations,
                   "Maximal iteration limit has to be larger than zero or -1 (unlimited).");
    PRECICE_ASSERT((maxIterations == INFINITE_MAX_ITERATIONS) || (minIterations <= maxIterations),
                   "Minimal iteration limit has to be smaller equal compared to the maximal iteration limit.");
  }

  if (math::equals(maxTime, UNDEFINED_MAX_TIME)) {
    _time = impl::TimeHandler();
  } else {
    _time = impl::TimeHandler(maxTime);
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

  profiling::Event e("waitAndSendData", profiling::Fundamental);

  for (const auto &data : sendData | boost::adaptors::map_values) {
    const auto &stamples = data->stamples();
    PRECICE_ASSERT(!stamples.empty());

    int nTimeSteps = data->timeStepsStorage().nTimes();
    PRECICE_ASSERT(nTimeSteps > 0);

    if (data->exchangeSubsteps()) {
      const Eigen::VectorXd timesAscending = data->timeStepsStorage().getTimes();
      sendNumberOfTimeSteps(m2n, nTimeSteps);
      sendTimes(m2n, timesAscending);

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
  profiling::Event e("waitAndReceiveData", profiling::Fundamental);
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
      time::Sample recvSample(data->getDimensions(), data->nVertices(), data->meshDimensions());
      m2n->receive(recvSample.values, data->getMeshID(), data->getDimensions());

      if (data->hasGradient()) {
        PRECICE_ASSERT(data->hasGradient());
        m2n->receive(recvSample.gradients, data->getMeshID(), data->getDimensions() * data->meshDimensions());
      }

      data->setSampleAtTime(getTime(), recvSample);
    }
  }
}

void BaseCouplingScheme::receiveDataForWindowEnd(const m2n::PtrM2N &m2n, const DataMap &receiveData)
{
  // @TODO This is a hack required until https://github.com/precice/precice/issues/1957
  // buffer current time
  const impl::TimeHandler oldTime = _time;
  // set _time state to point to end of this window such that _time in receiveData is at end of window
  _time.progressBy(_nextTimeWindowSize);
  // receive data for end of window
  this->receiveData(m2n, receiveData);
  // reset time state;
  _time = oldTime;
}

void BaseCouplingScheme::initializeWithZeroInitialData(const DataMap &receiveData)
{
  for (const auto &data : receiveData | boost::adaptors::map_values) {
    PRECICE_DEBUG("Initialize {} as zero.", data->getDataName());
    // just store already initialized zero sample to storage.
    data->setSampleAtTime(getTime(), data->sample());
  }
}

PtrCouplingData BaseCouplingScheme::addCouplingData(const mesh::PtrData &data, mesh::PtrMesh mesh, bool requiresInitialization, bool communicateSubsteps, CouplingData::Direction direction)
{
  const DataID id = data->getID();

  // new data
  if (_allData.count(id) == 0) {
    auto ptr = std::make_shared<CouplingData>(data, std::move(mesh), requiresInitialization, communicateSubsteps, direction);
    _allData.emplace(id, ptr);
    return ptr;
  }

  // data is already used by another exchange of this coupling scheme, use existing CouplingData
  auto &existing = _allData[id];
  PRECICE_CHECK(existing->getDirection() == direction,
                "Data \"{0}\" cannot be added for sending and for receiving. Please remove either <exchange data=\"{0}\" ... /> tag", data->getName());
  PRECICE_CHECK(direction == CouplingData::Direction::Send,
                "Data \"{0}\" cannot be received from multiple sources. Please remove either <exchange data=\"{0}\" ... /> tag", data->getName());
  PRECICE_CHECK(existing->exchangeSubsteps() == communicateSubsteps,
                "Data \"{0}\" is configured with substeps enabled and disabled at the same time. Please make the configuration consistent in the <exchange data=\"{0}\" ... substeps=\"True/False\" ... /> tags", data->getName());
  PRECICE_CHECK(existing->requiresInitialization == requiresInitialization,
                "Data \"{0}\" is configured with data initialization enabled and disabled at the same time. Please make the configuration consistent in the <exchange data=\"{0}\" ... /> tags", data->getName());

  return existing;
}

bool BaseCouplingScheme::isExplicitCouplingScheme() const
{
  PRECICE_ASSERT(_couplingMode != Undefined);
  return _couplingMode == Explicit;
}

void BaseCouplingScheme::setTimeWindowSize(double timeWindowSize)
{
  _timeWindowSize = timeWindowSize;
}

void BaseCouplingScheme::setNextTimeWindowSize(double timeWindowSize)
{
  _nextTimeWindowSize = timeWindowSize;
}

void BaseCouplingScheme::finalize()
{
  PRECICE_TRACE();
  checkCompletenessRequiredActions();
  PRECICE_ASSERT(_isInitialized, "Called finalize() before initialize().");
}

void BaseCouplingScheme::initialize()
{
  // initialize with zero data here, might eventually be overwritten in exchangeInitialData
  initializeReceiveDataStorage();
  // Initialize uses the template method pattern (https://en.wikipedia.org/wiki/Template_method_pattern).
  PRECICE_TRACE();
  PRECICE_ASSERT(not isInitialized());
  _time.resetTo(0);
  _timeWindows         = 1;
  _hasDataBeenReceived = false;

  checkCouplingDataAvailable();

  if (isImplicitCouplingScheme()) {
    storeIteration();
    if (not doesFirstStep()) {
      // reserve memory and initialize data with zero
      if (_acceleration) {
        _acceleration->initialize(getAccelerationData());
        if (auto qnAcceleration = std::dynamic_pointer_cast<precice::acceleration::BaseQNAcceleration>(_acceleration);
            qnAcceleration) {
          PRECICE_WARN_IF((qnAcceleration->getMaxUsedTimeWindows() == 0) && (qnAcceleration->getMaxUsedIterations() > _maxIterations), "The maximum number of iterations used in the quasi-Newton acceleration scheme is greater than the maximum number of iterations allowed in one time window. When time-windows-reused is set to 0, and therefore no previous time windows are reused, the actual max-used-ietrations is equal to max-iterations.");
        }
      }
    }
    requireAction(CouplingScheme::Action::WriteCheckpoint);
    initializeTXTWriters();
  }

  exchangeInitialData();

  _isInitialized = true;
}

void BaseCouplingScheme::reinitialize()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(isInitialized());

  if (isImplicitCouplingScheme()) {
    // overwrite past iteration with new samples
    for (const auto &data : _allData | boost::adaptors::map_values) {
      // TODO: reset CouplingData of changed meshes only #2102
      data->reinitialize();
    }

    if (not doesFirstStep()) {
      if (_acceleration) {
        _acceleration->initialize(getAccelerationData());
      }
    }
    initializeTXTWriters();
  }
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
        PRECICE_ASSERT(math::greater(getTime(), getWindowStartTime()));
        _timeWindows -= 1;
        _isTimeWindowComplete = false;
      } else { // write output, prepare for next window
        PRECICE_DEBUG("Convergence achieved");
        advanceTXTWriters();
        PRECICE_INFO("Time window completed");
        _isTimeWindowComplete = true;
        if (isCouplingOngoing()) {
          PRECICE_DEBUG("Setting require create checkpoint");
          requireAction(CouplingScheme::Action::WriteCheckpoint);
        }
      }
      // update iterations
      _totalIterations++;
      if (not hasConverged()) {
        _iterations++;
      } else {
        _iterations = 1;
      }
    } else {
      PRECICE_ASSERT(isExplicitCouplingScheme());
      PRECICE_INFO("Time window completed");
      _isTimeWindowComplete = true;
    }
    if (isCouplingOngoing()) {
      PRECICE_ASSERT(_hasDataBeenReceived);
    }

    // Update internal time tracking
    if (_isTimeWindowComplete) {
      _time.completeTimeWindow(_timeWindowSize);
    } else {
      _time.resetProgress();
    }
    _timeWindowSize = _nextTimeWindowSize;
  }
}

void BaseCouplingScheme::moveToNextWindow()
{
  PRECICE_TRACE(_timeWindows);
  for (auto &data : _allData | boost::adaptors::map_values) {
    data->moveToNextWindow();
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

double BaseCouplingScheme::getNextTimeWindowSize() const
{
  return _nextTimeWindowSize;
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

  // Check validness
  bool valid = math::greaterEquals(getNextTimeStepMaxSize(), timeToAdd);
  PRECICE_CHECK(valid,
                "The time step size given to preCICE in \"advance\" {} exceeds the maximum allowed time step size {} "
                "in the remaining of this time window. "
                "Did you restrict your time step size, \"dt = min(preciceDt, solverDt)\"? "
                "For more information, consult the adapter example in the preCICE documentation.",
                timeToAdd, getNextTimeStepMaxSize());

  // add time interval that has been computed in the solver to get the correct time remainder
  _time.progressBy(timeToAdd);

  return reachedEndOfTimeWindow();
}

bool BaseCouplingScheme::willDataBeExchanged(
    double lastSolverTimeStepSize) const
{
  PRECICE_TRACE(lastSolverTimeStepSize);
  double remainder = getNextTimeStepMaxSize() - lastSolverTimeStepSize;
  return not math::greater(remainder, 0.0);
}

bool BaseCouplingScheme::hasDataBeenReceived() const
{
  return _hasDataBeenReceived;
}

void BaseCouplingScheme::setDoesFirstStep(bool doesFirstStep)
{
  _doesFirstStep = doesFirstStep;
}

void BaseCouplingScheme::notifyDataHasBeenReceived()
{
  PRECICE_ASSERT(not _hasDataBeenReceived, "notifyDataHasBeenReceived() may only be called once within one coupling iteration. If this assertion is triggered this probably means that your coupling scheme has a bug.");
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
  return _time.time();
}

double BaseCouplingScheme::getTimeWindowStart() const
{
  return _time.windowStart();
}

int BaseCouplingScheme::getTimeWindows() const
{
  return _timeWindows;
}

double BaseCouplingScheme::getNextTimeStepMaxSize() const
{
  if (!isCouplingOngoing()) {
    return 0.0; // if coupling is not ongoing (i.e. coupling scheme reached end of window) the maximum time step size is zero
  }

  if (hasTimeWindowSize()) {
    return _time.untilWindowEnd(_timeWindowSize);
  } else {
    return _time.untilEnd();
  }
}

bool BaseCouplingScheme::isCouplingOngoing() const
{
  bool timestepsLeft = (_maxTimeWindows >= _timeWindows) || (_maxTimeWindows == UNDEFINED_TIME_WINDOWS);
  return timestepsLeft && !_time.reachedEnd();
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
  if (isCouplingOngoing()) {
    os << "iteration: " << _iterations; //_iterations;
    if ((_maxIterations != UNDEFINED_MAX_ITERATIONS) && (_maxIterations != INFINITE_MAX_ITERATIONS)) {
      os << " of " << _maxIterations;
    }
    if (_minIterations != UNDEFINED_MIN_ITERATIONS) {
      os << " (min " << _minIterations << ")";
    }
    os << ", ";
  }
  os << printBasicState(_timeWindows, getTime());
  std::string actionsState = printActionsState();
  if (!actionsState.empty()) {
    os << ", " << actionsState;
  }
  return os.str();
}

std::string BaseCouplingScheme::printBasicState(
    int    timeWindows,
    double time) const
{
  std::ostringstream os;
  if (isCouplingOngoing()) {
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
  } else {
    os << "Reached end at: final time-window: " << (timeWindows - 1) << ", final time: " << time;
  }
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

void BaseCouplingScheme::checkCouplingDataAvailable()
{
  PRECICE_TRACE();
  for (const auto &data : _allData | boost::adaptors::map_values) {
    if (data->getDirection() == CouplingData::Direction::Receive) {
      PRECICE_ASSERT(!data->stamples().empty(), "initializeReceiveDataStorage() didn't initialize data correctly");
    } else {
      PRECICE_CHECK(!data->stamples().empty(),
                    "Data {0} on mesh {1} doesn't contain any samples while initializing a coupling scheme of participant {2}. "
                    "There are two common configuration issues that may cause this. "
                    "Either, make sure participant {2} specifies data {0} to be written using tag <write-data mesh=\"{1}\" data=\"{0}\"/>. "
                    "Or ensure participant {2} defines a mapping to mesh {1} from a mesh using data {0}.",
                    data->getDataName(), data->getMeshName(), localParticipant());
    }
  }
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
    impl::PtrConvergenceMeasure measure)
{
  ConvergenceMeasureContext convMeasure;
  PRECICE_ASSERT(_allData.count(dataID) == 1, "Data with given data ID must exist!");
  convMeasure.couplingData = _allData.at(dataID);
  convMeasure.suffices     = suffices;
  convMeasure.strict       = strict;
  convMeasure.measure      = std::move(measure);
  _convergenceMeasures.push_back(convMeasure);
}

bool BaseCouplingScheme::measureConvergence()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not doesFirstStep());
  if (not utils::IntraComm::isSecondary()) {
    _convergenceWriter->writeData("TimeWindow", _timeWindows - 1);
    _convergenceWriter->writeData("Iteration", _iterations);
  }

  // If no convergence measures are defined, we never converge
  if (_convergenceMeasures.empty()) {
    PRECICE_INFO("No converge measures defined.");
    return false;
  }

  // There are convergence measures defined, so we need to check them
  bool allConverged = true;
  bool oneSuffices  = false; // at least one convergence measure suffices and did converge
  bool oneStrict    = false; // at least one convergence measure is strict and did not converge

  const bool reachedMinIterations = _iterations >= _minIterations;
  for (const auto &convMeasure : _convergenceMeasures) {
    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    PRECICE_ASSERT(convMeasure.measure.get() != nullptr);
    PRECICE_ASSERT(convMeasure.couplingData->previousIteration().size() == convMeasure.couplingData->values().size(), convMeasure.couplingData->previousIteration().size(), convMeasure.couplingData->values().size(), convMeasure.couplingData->getDataName());
    convMeasure.measure->measure(convMeasure.couplingData->previousIteration(), convMeasure.couplingData->values());

    if (not utils::IntraComm::isSecondary()) {
      _convergenceWriter->writeData(convMeasure.logHeader(), convMeasure.measure->getNormResidual());
    }

    if (not convMeasure.measure->isConvergence()) {
      allConverged = false;
      if (convMeasure.strict) {
        PRECICE_ASSERT(_maxIterations > 0);
        oneStrict = true;
        PRECICE_CHECK(_iterations < _maxIterations,
                      "The strict convergence measure for data \"" + convMeasure.couplingData->getDataName() +
                          "\" did not converge within the maximum allowed iterations, which terminates the simulation. "
                          "To avoid this forced termination do not mark the convergence measure as strict.");
      }
    } else if (convMeasure.suffices == true) {
      oneSuffices = true;
    }

    PRECICE_INFO(convMeasure.measure->printState(convMeasure.couplingData->getDataName()));
  }

  std::string messageSuffix;
  if (not reachedMinIterations) {
    messageSuffix = " but hasn't yet reached minimal amount of iterations";
  }
  if (allConverged) {
    PRECICE_INFO("All converged{}", messageSuffix);
  } else if (oneSuffices && not oneStrict) { // strict overrules suffices
    PRECICE_INFO("Sufficient measures converged{}", messageSuffix);
  }

  return reachedMinIterations && (allConverged || (oneSuffices && not oneStrict));
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
        _convergenceWriter->addData(convMeasure.logHeader(), io::TXTTableWriter::DOUBLE);
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
    bool converged = _iterations >= _minIterations && (_maxIterations < 0 || (_iterations < _maxIterations));
    _iterationsWriter->writeData("Convergence", converged ? 1 : 0);

    if (not doesFirstStep() && _acceleration) {
      std::shared_ptr<precice::acceleration::BaseQNAcceleration> qnAcceleration = std::dynamic_pointer_cast<precice::acceleration::BaseQNAcceleration>(_acceleration);
      if (qnAcceleration) {
        // Only write values for additional columns, if using a QN-based acceleration scheme
        _iterationsWriter->writeData("QNColumns", qnAcceleration->getLSSystemCols());
        _iterationsWriter->writeData("DeletedQNColumns", qnAcceleration->getDeletedColumns());
        _iterationsWriter->writeData("DroppedQNColumns", qnAcceleration->getDroppedColumns());
      } else {
        // non-breaking implementation uses "0" for these columns (delete columns in the future?)
        _iterationsWriter->writeData("QNColumns", 0);
        _iterationsWriter->writeData("DeletedQNColumns", 0);
        _iterationsWriter->writeData("DroppedQNColumns", 0);
      }
    }
  }
}

bool BaseCouplingScheme::reachedEndOfTimeWindow() const
{
  if (!hasTimeWindowSize()) {
    return true; // This participant will always do exactly one step to dictate second's time-window-size
  }

  return _time.reachedEndOfWindow(_timeWindowSize);
}

void BaseCouplingScheme::storeIteration()
{
  PRECICE_ASSERT(isImplicitCouplingScheme());

  for (const auto &data : _allData | boost::adaptors::map_values) {
    PRECICE_CHECK(!data->stamples().empty(),
                  "Data {0} on mesh {1} didn't contain any samples while attempting to send it to the coupling partner. "
                  "Make sure participant {2} specifies data {0} to be written using tag <write-data mesh=\"{1}\" data=\"{0}\"/>. "
                  "Alternatively, ensure participant {2} defines a mapping to mesh {1} from a mesh using data {0}.",
                  data->getDataName(), data->getMeshName(), localParticipant());
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
      profiling::Event e("accelerate", profiling::Fundamental);
      _acceleration->iterationsConverged(getAccelerationData(), getTimeWindowStart());
    }
    newConvergenceMeasurements();
  } else {
    // no convergence achieved for the coupling iteration within the current time window
    if (_acceleration) {
      profiling::Event e("accelerate", profiling::Fundamental);
      // Acceleration works on CouplingData::values(), so we retrieve the data from the storage, perform the acceleration and then put the data back into the storage. See also https://github.com/precice/precice/issues/1645.
      // @todo For acceleration schemes as described in "Rüth, B, Uekermann, B, Mehl, M, Birken, P, Monge, A, Bungartz, H-J. Quasi-Newton waveform iteration for partitioned surface-coupled multiphysics applications. https://doi.org/10.1002/nme.6443" we need a more elaborate implementation.

      // Load from storage into buffer
      for (auto &data : getAccelerationData() | boost::adaptors::map_values) {
        const auto &stamples = data->stamples();
        PRECICE_ASSERT(stamples.size() > 0);
        data->sample() = stamples.back().sample;
      }

      _acceleration->performAcceleration(getAccelerationData(), getTimeWindowStart());

      // Store from buffer
      // @todo Currently only data at end of window is accelerated. Remaining data in storage stays as it is.
      for (auto &data : getAccelerationData() | boost::adaptors::map_values) {
        data->setSampleAtTime(getTime(), data->sample());
      }
    }
  }
}

void BaseCouplingScheme::sendConvergence(const m2n::PtrM2N &m2n)
{
  PRECICE_ASSERT(isImplicitCouplingScheme());
  PRECICE_ASSERT(not doesFirstStep(), "For convergence information the sending participant is never the first one.");
  profiling::Event e("sendConvergence", profiling::Fundamental);
  m2n->send(_hasConverged);
}

void BaseCouplingScheme::receiveConvergence(const m2n::PtrM2N &m2n)
{
  PRECICE_ASSERT(isImplicitCouplingScheme());
  PRECICE_ASSERT(doesFirstStep(), "For convergence information the receiving participant is always the first one.");
  profiling::Event e("receiveConvergence", profiling::Fundamental);
  m2n->receive(_hasConverged);
}

double BaseCouplingScheme::getWindowStartTime() const
{
  return _time.windowStart();
}

double BaseCouplingScheme::getWindowEndTime() const
{
  return getWindowStartTime() + getTimeWindowSize();
}

bool BaseCouplingScheme::requiresSubsteps() const
{
  // Global toggle if a single send data uses substeps
  for (auto cpldata : _allData | boost::adaptors::map_values) {
    if (cpldata->getDirection() == CouplingData::Direction::Send && cpldata->exchangeSubsteps()) {
      return true;
    }
  }
  return false;
}

ImplicitData BaseCouplingScheme::implicitDataToReceive() const
{
  // This implementation covers everything except SerialImplicit
  if (!isImplicitCouplingScheme()) {
    return {};
  }

  ImplicitData idata;
  for (auto cpldata : _allData | boost::adaptors::map_values) {
    if (cpldata->getDirection() == CouplingData::Direction::Receive) {
      idata.add(cpldata->getDataID(), false);
    }
  }
  return idata;
}

std::string BaseCouplingScheme::localParticipant() const
{
  return _localParticipant;
}

} // namespace precice::cplscheme
