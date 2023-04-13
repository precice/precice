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
  PRECICE_ASSERT(m2n.get() != nullptr);
  PRECICE_ASSERT(m2n->isConnected());

  for (const auto &data : sendData | boost::adaptors::map_values) {
    auto timesAscending = data->getStoredTimesAscending();

    PRECICE_ASSERT(math::equals(timesAscending(0), time::Storage::WINDOW_START), timesAscending(0));                                               // assert that first element is time::Storage::WINDOW_START
    PRECICE_ASSERT(math::equals(timesAscending(timesAscending.size() - 1), time::Storage::WINDOW_END), timesAscending(timesAscending.size() - 1)); // assert that last element is time::Storage::WINDOW_END

    sendNumberOfTimeSteps(m2n, timesAscending.size());
    sendTimes(m2n, timesAscending);

    auto serializedSamples = data->getSerialized();
    // Data is actually only send if size>0, which is checked in the derived classes implementation
    m2n->send(serializedSamples, data->getMeshID(), data->getDimensions() * timesAscending.size());

    if (data->hasGradient()) {
      m2n->send(data->gradientValues(), data->getMeshID(), data->getDimensions() * data->meshDimensions());
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
  auto times = Eigen::VectorXd(nTimeSteps);
  for (int i = 0; i < nTimeSteps; i++) {
    m2n->receive(times(i));
  }
  PRECICE_DEBUG("Received times {}", times);
  return times;
}

void BaseCouplingScheme::receiveData(const m2n::PtrM2N &m2n, const DataMap &receiveData)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(m2n.get());
  PRECICE_ASSERT(m2n->isConnected());
  for (const auto &data : receiveData | boost::adaptors::map_values) {
    PRECICE_DEBUG("BaseCouplingScheme::receiveData of DataID {}", data->getDataID());

    int nTimeSteps = receiveNumberOfTimeSteps(m2n);

    PRECICE_ASSERT(nTimeSteps > 0);
    auto timesAscending = receiveTimes(m2n, nTimeSteps);

    PRECICE_ASSERT(math::equals(timesAscending(0), time::Storage::WINDOW_START), timesAscending(0));                                               // assert that first element is time::Storage::WINDOW_START
    PRECICE_ASSERT(math::equals(timesAscending(timesAscending.size() - 1), time::Storage::WINDOW_END), timesAscending(timesAscending.size() - 1)); // assert that last element is time::Storage::WINDOW_END

    auto serializedSamples = Eigen::VectorXd(nTimeSteps * data->getSize());
    // Data is only received on ranks with size>0, which is checked in the derived class implementation
    m2n->receive(serializedSamples, data->getMeshID(), data->getDimensions() * nTimeSteps);

    data->storeFromSerialized(timesAscending, serializedSamples);

    if (data->hasGradient()) {
      m2n->receive(data->gradientValues(), data->getMeshID(), data->getDimensions() * data->meshDimensions());
    }
  }
}

void BaseCouplingScheme::initializeZeroReceiveData(const DataMap &receiveData)
{
  for (const auto &data : receiveData | boost::adaptors::map_values) {
    auto zeroData = Eigen::VectorXd::Zero(data->getSize());
    data->initializeStorage(zeroData);
  }
}

PtrCouplingData BaseCouplingScheme::addCouplingData(const mesh::PtrData &data, mesh::PtrMesh mesh, bool requiresInitialization)
{
  int             id = data->getID();
  PtrCouplingData ptrCplData;
  if (!utils::contained(id, _allData)) { // data is not used by this coupling scheme yet, create new CouplingData
    if (isExplicitCouplingScheme()) {
      /// ptrCplData = std::make_shared<CouplingData>(data, std::move(mesh), requiresInitialization, precice::cplscheme::CouplingSchemeConfiguration::getWaveformUsedOrder(this->_localParticipant, data->getName()));
      ptrCplData = std::make_shared<CouplingData>(data, std::move(mesh), requiresInitialization, 1);
    } else {
      ptrCplData = std::make_shared<CouplingData>(data, std::move(mesh), requiresInitialization, 1, getExtrapolationOrder());
    }
    _allData.emplace(id, ptrCplData);
  } else { // data is already used by another exchange of this coupling scheme, use existing CouplingData
    ptrCplData = _allData[id];
  }
  return ptrCplData;
}

bool BaseCouplingScheme::isExplicitCouplingScheme()
{
  PRECICE_ASSERT(_couplingMode != Undefined);
  return _couplingMode == Explicit;
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

  //Need to initialize the waveforms before I store them.
  initializeSendDataStorage();
  exchangeInitialData();

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
      storeSendValuesAtTime(relativeDt);
    } else {
      // We don't support subcycling here, because this is complicated. Therefore, use same strategy like for explicit coupling and just use a single value at end of window.
      // Possible solution: Don't scale times to [0,1], but leave them as they are. Then we would also allow times > 1. We then have two options:
      // 1) scale the times back later when the time window size is known (to still benefit from the simpler handling, if all times are scaled to [0,1]).
      // 2) generally use times in the interval [0, timeWindowSize]. This makes the implementation probably a bit more complicated, but also more consistent.
      if (reachedEndOfTimeWindow()) {                     // only necessary to trigger at end of time window.
        storeSendValuesAtTime(time::Storage::WINDOW_END); // only write data at end of window
      }
    }
  } else {
    // work-around for explicit coupling, because it does not support waveform relaxation.
    if (reachedEndOfTimeWindow()) {                     // only necessary to trigger at end of time window.
      storeSendValuesAtTime(time::Storage::WINDOW_END); // only write data at end of window
    }
  }

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
  PRECICE_TRACE(_timeWindows, _time);
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
        _time = _time - _computedTimeWindowPart;
        _timeWindows -= 1;
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

bool BaseCouplingScheme::isInitialized() const
{
  return _isInitialized;
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

double BaseCouplingScheme::getComputedTimeWindowPart()
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
  return math::equals(getThisTimeWindowRemainder(), 0.0, _eps);
}

void BaseCouplingScheme::storeIteration()
{
  PRECICE_ASSERT(isImplicitCouplingScheme());
  for (const auto &data : _allData | boost::adaptors::map_values) {
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

int BaseCouplingScheme::getExtrapolationOrder()
{
  return _extrapolationOrder;
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
      _acceleration->iterationsConverged(getAccelerationData());
    }
    newConvergenceMeasurements();
  } else {
    // no convergence achieved for the coupling iteration within the current time window
    if (_acceleration) {
      /**
       * Acceleration works on CouplingData::values(), so we retreive the data from the storage, perform the acceleration and then put the data back into the storage.
       *
       * There are generally two possibilities:
       *
       * 1) Only overwrite the data in CouplingData::_timeStepsStorage that is part of the receive data
       * 2) overwrite the data in CouplingData::_timeStepsStorage for all data in getAccelerationData()
       *
       * We are using strategy 2), because it's easier to access getAccelerationData() from here. However, this also means that
       * we have to make sure that the send data is stored in the CouplingData::_timeStepsStorage - even though this is not
       * needed at the moment. Important note: In https://github.com/precice/precice/pull/1414 also send data requires to keep
       * track of _timeStepsStorage for subcycling. So it will become simpler as soon as subcycling is fully implemented.
       */
      // @todo For other Acceleration schemes as described in "Rüth, B, Uekermann, B, Mehl, M, Birken, P, Monge, A, Bungartz, H-J. Quasi-Newton waveform iteration for partitioned surface-coupled multiphysics applications. Int J Numer Methods Eng. 2021; 122: 5236– 5257. https://doi.org/10.1002/nme.6443" we need a more elaborate implementation.
      // Put values into data->values() for acceleration
      // @todo breaks for CplSchemeTests/ParallelImplicitCouplingSchemeTests. Why? @fsimonis
      // for (auto &data : getAccelerationData() | boost::adaptors::map_values) {
      //   data->values() = data->getValuesAtTime(time::Storage::WINDOW_END);
      // }
      for (auto &pair : getAccelerationData()) {
        pair.second->values() = pair.second->getValuesAtTime(time::Storage::WINDOW_END);
      }

      _acceleration->performAcceleration(getAccelerationData());
      // Need to copy over the accelerated data into values again
      for (auto &pair : getAccelerationData()) {
        pair.second->values() = pair.second->getValuesAtTime(time::Storage::WINDOW_END);
      }
      // @todo breaks for CplSchemeTests/ParallelImplicitCouplingSchemeTests. Why? @fsimonis
      // for (auto &data : getAccelerationData() | boost::adaptors::map_values) {
      //   data->storeValuesAtTime(time::Storage::WINDOW_END, data->values(), mustOverwrite);
      // }
      // The acceleration is applied to the data in the acceleration class not here.
      // for (auto &pair : getAccelerationData()) {
      //   bool mustOverwrite = true;
      //   pair.second->storeValuesAtTime(time::Storage::WINDOW_END, pair.second->values(), mustOverwrite); // @todo: Here might be an error. Do we actually apply any acceleration to the values that are not at WINDOW_END?
      // }
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
