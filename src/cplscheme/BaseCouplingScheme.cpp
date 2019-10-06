#include "BaseCouplingScheme.hpp"
#include <Eigen/Core>
#include <limits>
#include <sstream>
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "acceleration/Acceleration.hpp"
#include "io/TXTReader.hpp"
#include "io/TXTWriter.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "math/math.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"

namespace precice
{
namespace cplscheme
{

BaseCouplingScheme::BaseCouplingScheme(
    double maxTime,
    int    maxTimesteps,
    double timestepLength,
    int    validDigits)
    : _couplingMode(Undefined),
      _firstParticipant("unknown"),
      _secondParticipant("unknown"),
      _localParticipant("unknown"),
      _eps(std::pow(10.0, -1 * validDigits)),
      _iterationsCoarseOptimization(-1),
      _maxTime(maxTime),
      _maxTimesteps(maxTimesteps),
      _iterations(-1),
      _totalIterationsCoarseOptimization(-1),
      _maxIterations(-1),
      _totalIterations(-1),
      _timesteps(0),
      _timestepLength(timestepLength),
      _validDigits(validDigits)
{
  PRECICE_CHECK(not((maxTime != UNDEFINED_TIME) && (maxTime < 0.0)),
        "Maximum time has to be larger than zero!");
  PRECICE_CHECK(not((maxTimesteps != UNDEFINED_TIMESTEPS) && (maxTimesteps < 0)),
        "Maximum timestep number has to be larger than zero!");
  PRECICE_CHECK(not((timestepLength != UNDEFINED_TIMESTEP_LENGTH) && (timestepLength < 0.0)),
        "Timestep length has to be larger than zero!");
  PRECICE_CHECK((_validDigits >= 1) && (_validDigits < 17),
        "Valid digits of timestep length has to be between 1 and 16!");
}

BaseCouplingScheme::BaseCouplingScheme(
    double                        maxTime,
    int                           maxTimesteps,
    double                        timestepLength,
    int                           validDigits,
    const std::string &           firstParticipant,
    const std::string &           secondParticipant,
    const std::string &           localParticipant,
    m2n::PtrM2N                   m2n,
    int                           maxIterations,
    constants::TimesteppingMethod dtMethod)
  :   _firstParticipant(firstParticipant),
      _secondParticipant(secondParticipant),
      _localParticipant(localParticipant),
      _eps(std::pow(10.0, -1 * validDigits)),
      _iterationsCoarseOptimization(1),
      _m2n(m2n),
      _maxTime(maxTime),
      _maxTimesteps(maxTimesteps),
      _iterations(1),
      _totalIterationsCoarseOptimization(1),
      _maxIterations(maxIterations),
      _totalIterations(1),
      _timesteps(1),
      _timestepLength(timestepLength),
      _validDigits(validDigits)
{
  PRECICE_CHECK(not((maxTime != UNDEFINED_TIME) && (maxTime < 0.0)),
        "Maximum time has to be larger than zero!");
  PRECICE_CHECK(not((maxTimesteps != UNDEFINED_TIMESTEPS) && (maxTimesteps < 0)),
        "Maximum timestep number has to be larger than zero!");
  PRECICE_CHECK(not((timestepLength != UNDEFINED_TIMESTEP_LENGTH) && (timestepLength < 0.0)),
        "Timestep length has to be larger than zero!");
  PRECICE_CHECK((_validDigits >= 1) && (_validDigits < 17),
        "Valid digits of timestep length has to be between 1 and 16!");
  PRECICE_CHECK(_firstParticipant != _secondParticipant,
        "First participant and second participant must have different names! Called from BaseCoupling.");
  if (dtMethod == constants::FIXED_DT) {
    PRECICE_CHECK(hasTimestepLength(),
          "Timestep length value has to be given when the fixed timestep length method "
          << "is chosen for an implicit coupling scheme!");
  }
  if (localParticipant == _firstParticipant) {
    _doesFirstStep = true;
    if (dtMethod == constants::FIRST_PARTICIPANT_SETS_DT) {
      _participantSetsDt = true;
      setTimestepLength(UNDEFINED_TIMESTEP_LENGTH);
    }
  } else if (localParticipant == _secondParticipant) {
    if (dtMethod == constants::FIRST_PARTICIPANT_SETS_DT) {
      _participantReceivesDt = true;
    }
  } else {
    PRECICE_ERROR("Name of local participant \""
          << localParticipant << "\" does not match any "
          << "participant specified for the coupling scheme!");
  }
  PRECICE_CHECK((maxIterations > 0) || (maxIterations == -1),
        "Maximal iteration limit has to be larger than zero!");
}

void BaseCouplingScheme::receiveAndSetDt()
{
  PRECICE_TRACE();
  if (participantReceivesDt()) {
    double dt = UNDEFINED_TIMESTEP_LENGTH;
    getM2N()->receive(dt);
    PRECICE_DEBUG("Received timestep length of " << dt);
    PRECICE_ASSERT(not math::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
    setTimestepLength(dt);
  }
}

void BaseCouplingScheme::sendDt()
{
  PRECICE_TRACE();
  if (participantSetsDt()) {
    PRECICE_DEBUG("sending timestep length of " << getComputedTimestepPart());
    getM2N()->send(getComputedTimestepPart());
  }
}

void BaseCouplingScheme::addDataToSend(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize)
{
  PRECICE_TRACE();
  int id = data->getID();
  if (!utils::contained(id, _sendData)) {
    PtrCouplingData     ptrCplData(new CouplingData(&(data->values()), mesh, initialize, data->getDimensions()));
    DataMap::value_type pair = std::make_pair(id, ptrCplData);
    _sendData.insert(pair);
  } else {
    PRECICE_ERROR("Data \"" << data->getName() << "\" cannot be added twice for sending!");
  }
}

void BaseCouplingScheme::addDataToReceive(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize)
{
  PRECICE_TRACE();
  int id = data->getID();
  if (!utils::contained(id, _receiveData)) {
    PtrCouplingData     ptrCplData(new CouplingData(&(data->values()), mesh, initialize, data->getDimensions()));
    DataMap::value_type pair = std::make_pair(id, ptrCplData);
    _receiveData.insert(pair);
  } else {
    PRECICE_ERROR("Data \"" << data->getName() << "\" cannot be added twice for receiving!");
  }
}

void BaseCouplingScheme::sendState(
    com::PtrCommunication communication,
    int                   rankReceiver)
{
  PRECICE_TRACE(rankReceiver);
  PRECICE_ASSERT(communication.get() != nullptr);
  PRECICE_ASSERT(communication->isConnected());
  communication->send(_maxTime, rankReceiver);
  communication->send(_maxTimesteps, rankReceiver);
  communication->send(_timestepLength, rankReceiver);
  communication->send(_time, rankReceiver);
  communication->send(_timesteps, rankReceiver);
  communication->send(_computedTimestepPart, rankReceiver);
  //communication->send(_maxLengthNextTimestep, rankReceiver);
  communication->send(_isInitialized, rankReceiver);
  communication->send(_isCouplingTimestepComplete, rankReceiver);
  communication->send(_hasDataBeenExchanged, rankReceiver);
  communication->send((int) _actions.size(), rankReceiver);
  for (const std::string &action : _actions) {
    communication->send(action, rankReceiver);
  }
  communication->send(_maxIterations, rankReceiver);
  communication->send(_iterations, rankReceiver);
  communication->send(_iterationsCoarseOptimization, rankReceiver); // new, correct?? TODO
  communication->send(_totalIterations, rankReceiver);
}

void BaseCouplingScheme::receiveState(
    com::PtrCommunication communication,
    int                   rankSender)
{
  PRECICE_TRACE(rankSender);
  PRECICE_ASSERT(communication.get() != nullptr);
  PRECICE_ASSERT(communication->isConnected());
  communication->receive(_maxTime, rankSender);
  communication->receive(_maxTimesteps, rankSender);
  communication->receive(_timestepLength, rankSender);
  communication->receive(_time, rankSender);
  communication->receive(_timesteps, rankSender);
  communication->receive(_computedTimestepPart, rankSender);
  //communication->receive(_maxLengthNextTimestep, rankSender);
  communication->receive(_isInitialized, rankSender);
  communication->receive(_isCouplingTimestepComplete, rankSender);
  communication->receive(_hasDataBeenExchanged, rankSender);
  int actionsSize = 0;
  communication->receive(actionsSize, rankSender);
  _actions.clear();
  for (int i = 0; i < actionsSize; i++) {
    std::string action;
    communication->receive(action, rankSender);
    _actions.insert(action);
  }
  communication->receive(_maxIterations, rankSender);
  int subIteration = -1;
  communication->receive(subIteration, rankSender);
  _iterations = subIteration;
  communication->receive(subIteration, rankSender); // new, correct?? TODO
  _iterationsCoarseOptimization = subIteration;     // new, correct? TODO
  communication->receive(_totalIterations, rankSender);
}

std::vector<int> BaseCouplingScheme::sendData(m2n::PtrM2N m2n)
{
  PRECICE_TRACE();

  std::vector<int> sentDataIDs;
  PRECICE_ASSERT(m2n.get() != nullptr);
  PRECICE_ASSERT(m2n->isConnected());
  for (const DataMap::value_type &pair : _sendData) {
    //std::cout<<"\nsend data id="<<pair.first<<": "<<*(pair.second->values)<<'\n';
    int size = pair.second->values->size();
    m2n->send(pair.second->values->data(), size, pair.second->mesh->getID(), pair.second->dimension);
    sentDataIDs.push_back(pair.first);
  }
  PRECICE_DEBUG("Number of sent data sets = " << sentDataIDs.size());
  return sentDataIDs;
}

std::vector<int> BaseCouplingScheme::receiveData(
    m2n::PtrM2N m2n)
{
  PRECICE_TRACE();
  std::vector<int> receivedDataIDs;
  PRECICE_ASSERT(m2n.get() != nullptr);
  PRECICE_ASSERT(m2n->isConnected());

  for (DataMap::value_type &pair : _receiveData) {
    int size = pair.second->values->size();
    //std::cout<<"\nreceive data id="<<pair.first<<": "<<*(pair.second->values)<<'\n';
    m2n->receive(pair.second->values->data(), size, pair.second->mesh->getID(), pair.second->dimension);
    receivedDataIDs.push_back(pair.first);
  }
  PRECICE_DEBUG("Number of received data sets = " << receivedDataIDs.size());

  return receivedDataIDs;
}

int BaseCouplingScheme::getVertexOffset(
    std::map<int, int> &vertexDistribution,
    int                 rank,
    int                 dim)
{
  int sum = 0;
  for (int i = 0; i < rank; i++) {
    sum += vertexDistribution[i];
  }
  return sum * dim;
}

CouplingData *BaseCouplingScheme::getSendData(
    int dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _sendData.find(dataID);
  if (iter != _sendData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

CouplingData *BaseCouplingScheme::getReceiveData(
    int dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _receiveData.find(dataID);
  if (iter != _receiveData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

void BaseCouplingScheme::finalize()
{
  PRECICE_TRACE();
  checkCompletenessRequiredActions();
  PRECICE_CHECK(isInitialized(), "Called finalize() before initialize()!");
}

void BaseCouplingScheme::setExtrapolationOrder(
    int order)
{
  PRECICE_CHECK((order == 0) || (order == 1) || (order == 2),
        "Extrapolation order has to be  0, 1, or 2!");
  _extrapolationOrder = order;
}

// @todo extrapolation of data should only be done for the fine cplData -> then copied to the coarse cplData
void BaseCouplingScheme::extrapolateData(DataMap &data)
{
  PRECICE_TRACE(_timesteps);
  if ((_extrapolationOrder == 1) || getTimesteps() == 2) { //timesteps is increased before extrapolate is called
    PRECICE_INFO("Performing first order extrapolation");
    for (DataMap::value_type &pair : data) {
      PRECICE_DEBUG("Extrapolate data: " << pair.first);
      PRECICE_ASSERT(pair.second->oldValues.cols() > 1);
      Eigen::VectorXd &values       = *pair.second->values;
      pair.second->oldValues.col(0) = values;  // = x^t
      values *= 2.0;                           // = 2*x^t
      values -= pair.second->oldValues.col(1); // = 2*x^t - x^(t-1)
      utils::shiftSetFirst(pair.second->oldValues, values);
    }
  } else if (_extrapolationOrder == 2) {
    PRECICE_INFO("Performing second order extrapolation");
    for (DataMap::value_type &pair : data) {
      PRECICE_ASSERT(pair.second->oldValues.cols() > 2);
      Eigen::VectorXd &values     = *pair.second->values;
      auto             valuesOld1 = pair.second->oldValues.col(1);
      auto             valuesOld2 = pair.second->oldValues.col(2);

      pair.second->oldValues.col(0) = values; // = x^t
      values *= 2.5;                          // = 2.5 x^t
      values -= valuesOld1 * 2.0;             // = 2.5x^t - 2x^(t-1)
      values += valuesOld2 * 0.5;             // = 2.5x^t - 2x^(t-1) + 0.5x^(t-2)
      utils::shiftSetFirst(pair.second->oldValues, values);
    }
  } else {
    PRECICE_ERROR("Called extrapolation with order != 1,2!");
  }
}

bool BaseCouplingScheme::hasTimestepLength() const
{
  return not math::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH);
}

double BaseCouplingScheme::getTimestepLength() const
{
  PRECICE_ASSERT(not math::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH));
  return _timestepLength;
}

void BaseCouplingScheme::addComputedTime(
    double timeToAdd)
{
  PRECICE_TRACE(timeToAdd, _time);
  PRECICE_ASSERT(isCouplingOngoing(), "Invalid call of addComputedTime() after simulation end!");

  // add time interval that has been computed in the solver to get the correct time remainder
  _computedTimestepPart += timeToAdd;
  _time += timeToAdd;

  // Check validness
  bool valid = math::greaterEquals(getThisTimestepRemainder(), 0.0, _eps);
  PRECICE_CHECK(valid, "The computed timestep length of "
                   << timeToAdd << " exceeds the maximum timestep limit of "
                   << _timestepLength - _computedTimestepPart + timeToAdd
                   << " for this time step!");
}

bool BaseCouplingScheme::willDataBeExchanged(
    double lastSolverTimestepLength) const
{
  PRECICE_TRACE(lastSolverTimestepLength);
  double remainder = getThisTimestepRemainder() - lastSolverTimestepLength;
  return not math::greater(remainder, 0.0, _eps);
}

bool BaseCouplingScheme::hasDataBeenExchanged() const
{
  return _hasDataBeenExchanged;
}

void BaseCouplingScheme::setHasDataBeenExchanged(
    bool hasDataBeenExchanged)
{
  _hasDataBeenExchanged = hasDataBeenExchanged;
}

double BaseCouplingScheme::getTime() const
{
  return _time;
}

int BaseCouplingScheme::getTimesteps() const
{
  return _timesteps;
}

std::vector<std::string> BaseCouplingScheme::getCouplingPartners() const
{
  std::vector<std::string> partnerNames;
  // Add non-local participant
  if (doesFirstStep()) {
    partnerNames.push_back(_secondParticipant);
  } else {
    partnerNames.push_back(_firstParticipant);
  }
  return partnerNames;
}

double BaseCouplingScheme::getThisTimestepRemainder() const
{
  PRECICE_TRACE();
  double remainder = 0.0;
  if (not math::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH)) {
    remainder = _timestepLength - _computedTimestepPart;
  }
  PRECICE_DEBUG("return " << remainder);
  return remainder;
}

double BaseCouplingScheme::getNextTimestepMaxLength() const
{
  if (math::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH)) {
    if (math::equals(_maxTime, UNDEFINED_TIME)) {
      return std::numeric_limits<double>::max();
    } else {
      return _maxTime - _time;
    }
  }
  return _timestepLength - _computedTimestepPart;
}

bool BaseCouplingScheme::isCouplingOngoing() const
{
  bool timeLeft      = math::greater(_maxTime, _time, _eps) || math::equals(_maxTime, UNDEFINED_TIME);
  bool timestepsLeft = (_maxTimesteps >= _timesteps) || (_maxTimesteps == UNDEFINED_TIMESTEPS);
  return timeLeft && timestepsLeft;
}

bool BaseCouplingScheme::isCouplingTimestepComplete() const
{
  return _isCouplingTimestepComplete;
}

bool BaseCouplingScheme::isActionRequired(
    const std::string &actionName) const
{
  return _actions.count(actionName) > 0;
}

void BaseCouplingScheme::performedAction(
    const std::string &actionName)
{
  _actions.erase(actionName);
}

void BaseCouplingScheme::requireAction(
    const std::string &actionName)
{
  _actions.insert(actionName);
}

/// @todo: insert _iterationsCoarseOptimization in print state
std::string BaseCouplingScheme::printCouplingState() const
{
  std::ostringstream os;
  os << "it " << _iterations; //_iterations;
  if (getMaxIterations() != -1) {
    os << " of " << getMaxIterations();
  }
  os << " | " << printBasicState(_timesteps, _time) << " | " << printActionsState();
  return os.str();
}

std::string BaseCouplingScheme::printBasicState() const
{
  std::ostringstream os;
  os << printBasicState(_timesteps, _time);
  return os.str();
}

std::string BaseCouplingScheme::printBasicState(
    int    timesteps,
    double time) const
{
  std::ostringstream os;
  os << "dt# " << timesteps;
  if (_maxTimesteps != UNDEFINED_TIMESTEPS) {
    os << " of " << _maxTimesteps;
  }
  os << " | t " << time;
  if (_maxTime != UNDEFINED_TIME) {
    os << " of " << _maxTime;
  }
  if (_timestepLength != UNDEFINED_TIMESTEP_LENGTH) {
    os << " | dt " << _timestepLength;
  }
  if ((_timestepLength != UNDEFINED_TIMESTEP_LENGTH) || (_maxTime != UNDEFINED_TIME)) {
    os << " | max dt " << getNextTimestepMaxLength();
  }
  os << " | ongoing ";
  isCouplingOngoing() ? os << "yes" : os << "no";
  os << " | dt complete ";
  _isCouplingTimestepComplete ? os << "yes" : os << "no";
  return os.str();
}

std::string BaseCouplingScheme::printActionsState() const
{
  std::ostringstream os;
  for (const std::string &actionName : _actions) {
    os << actionName << " | ";
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
    PRECICE_ERROR("Unfulfilled required actions: " << stream.str() << "!");
  }
}

int BaseCouplingScheme::getValidDigits() const
{
  return _validDigits;
}

void BaseCouplingScheme::setupDataMatrices(DataMap &data)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Data size: " << data.size());
  // Reserve storage for convergence measurement of send and receive data values
  for (ConvergenceMeasure &convMeasure : _convergenceMeasures) {
    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    if (convMeasure.couplingData->oldValues.cols() < 1) {
      utils::append(convMeasure.couplingData->oldValues,
                    (Eigen::MatrixXd) Eigen::MatrixXd::Zero(convMeasure.couplingData->values->size(), 1));
    }
  }
  // Reserve storage for extrapolation of data values
  if (_extrapolationOrder > 0) {
    for (DataMap::value_type &pair : data) {
      int cols = pair.second->oldValues.cols();
      PRECICE_DEBUG("Add cols: " << pair.first << ", cols: " << cols);
      PRECICE_ASSERT(cols <= 1, cols);
      utils::append(pair.second->oldValues,
                    (Eigen::MatrixXd) Eigen::MatrixXd::Zero(pair.second->values->size(), _extrapolationOrder + 1 - cols));
    }
  }
}

void BaseCouplingScheme::setIterationAcceleration(
    acceleration::PtrAcceleration acceleration)
{
  PRECICE_ASSERT(acceleration.get() != nullptr);
  _acceleration = acceleration;

  // if multilevel based approach, i.e., manifold mapping, we have to start
  // with the evaluation/optimization of the coarse model representation.
  // otherwise, we start with the fine model representation as it's the only one
  _isCoarseModelOptimizationActive = _acceleration->isMultilevelBasedApproach();
  if (_acceleration->isMultilevelBasedApproach()) {
    _acceleration->setCoarseModelOptimizationActive(&_isCoarseModelOptimizationActive);
    // also initialize the iteration counters with 0, as scheme starts with coarse model evaluation
    _iterations      = 0;
    _totalIterations = 0;
  }
}

void BaseCouplingScheme::setupConvergenceMeasures()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not doesFirstStep());
  PRECICE_CHECK(not _convergenceMeasures.empty(),
        "At least one convergence measure has to be defined for "
            << "an implicit coupling scheme!");
  for (ConvergenceMeasure &convMeasure : _convergenceMeasures) {
    int dataID = convMeasure.data->getID();
    if ((getSendData(dataID) != nullptr)) {
      convMeasure.couplingData = getSendData(dataID);
    } else {
      convMeasure.couplingData = getReceiveData(dataID);
      PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    }
  }
}

void BaseCouplingScheme::newConvergenceMeasurements()
{
  PRECICE_TRACE();
  for (ConvergenceMeasure &convMeasure : _convergenceMeasures) {
    PRECICE_ASSERT(convMeasure.measure.get() != nullptr);
    convMeasure.measure->newMeasurementSeries();
  }
}

void BaseCouplingScheme::addConvergenceMeasure(
    mesh::PtrData               data,
    bool                        suffices,
    int                         level,
    impl::PtrConvergenceMeasure measure)
{
  ConvergenceMeasure convMeasure;
  convMeasure.data     = std::move(data);
  convMeasure.couplingData = nullptr;
  convMeasure.suffices = suffices;
  convMeasure.level    = level;
  convMeasure.measure  = std::move(measure);
  _convergenceMeasures.push_back(convMeasure);
  _firstResiduumNorm.push_back(0);
}

bool BaseCouplingScheme::measureConvergence(
    std::map<int, Eigen::VectorXd> &designSpecifications)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not doesFirstStep());
  bool allConverged = true;
  bool oneSuffices  = false;
  PRECICE_ASSERT(_convergenceMeasures.size() > 0);
  if (not utils::MasterSlave::isSlave()) {
    _convergenceWriter->writeData("Timestep", _timesteps);
    _convergenceWriter->writeData("Iteration", _iterations);
  }
  for (size_t i = 0; i < _convergenceMeasures.size(); i++) {
    ConvergenceMeasure &convMeasure = _convergenceMeasures[i];

    // only apply convergence measures for fine model optimization, i.e., coupling
    if (convMeasure.level > 0)
      continue;

    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    PRECICE_ASSERT(convMeasure.measure.get() != nullptr);
    const auto &    oldValues = convMeasure.couplingData->oldValues.col(0);
    Eigen::VectorXd q         = Eigen::VectorXd::Zero(convMeasure.couplingData->values->size());
    if (designSpecifications.find(convMeasure.data->getID()) != designSpecifications.end())
      q = designSpecifications.at(convMeasure.data->getID());

    convMeasure.measure->measure(oldValues, *convMeasure.couplingData->values, q);

    if (not utils::MasterSlave::isSlave()) {
      std::stringstream sstm;
      sstm << "resNorm(" << i << ")";
      _convergenceWriter->writeData(sstm.str(), convMeasure.measure->getNormResidual());
    }

    if (_iterations == 1)
      _firstResiduumNorm[i] = convMeasure.measure->getNormResidual();

    if (not convMeasure.measure->isConvergence()) {
      allConverged = false;
    } else if (convMeasure.suffices == true) {
      oneSuffices = true;
    }
    PRECICE_INFO(convMeasure.measure->printState());
  }

  if (allConverged) {
    PRECICE_INFO("All converged");
  } else if (oneSuffices) {
    PRECICE_INFO("Sufficient measure converged");
  }

  return allConverged || oneSuffices;
}

/// @todo: ugly hack with design specifications, however, getting them here is not possible as
// parallel coupling scheme and multi-coupling scheme  need allData and not only getSendData()
bool BaseCouplingScheme::measureConvergenceCoarseModelOptimization(
    std::map<int, Eigen::VectorXd> &designSpecifications)
{
  PRECICE_TRACE();
  bool allConverged = true;
  bool oneSuffices  = false;
  PRECICE_ASSERT(_convergenceMeasures.size() > 0);
  for (ConvergenceMeasure &convMeasure : _convergenceMeasures) {

    // only apply convergence measures for coarse model optimization
    if (convMeasure.level == 0)
      continue;

    std::cout << "  measure convergence coarse measure, data:" << convMeasure.data->getName() << '\n';
    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
    PRECICE_ASSERT(convMeasure.measure.get() != nullptr);
    const auto &    oldValues = convMeasure.couplingData->oldValues.col(0);
    Eigen::VectorXd q         = Eigen::VectorXd::Zero(convMeasure.couplingData->values->size());
    if (designSpecifications.find(convMeasure.data->getID()) != designSpecifications.end())
      q = designSpecifications.at(convMeasure.data->getID());

    convMeasure.measure->measure(oldValues, *convMeasure.couplingData->values, q);

    if (not convMeasure.measure->isConvergence()) {
      allConverged = false;
    } else if (convMeasure.suffices == true) {
      oneSuffices = true;
    }
    PRECICE_INFO('<' << convMeasure.data->getName() << '<' << convMeasure.measure->printState());
  }

  if (allConverged) {
    PRECICE_INFO("All converged");
  } else if (oneSuffices) {
    PRECICE_INFO("Sufficient measure converged");
  }

  return allConverged || oneSuffices;
}

void BaseCouplingScheme::initializeTXTWriters()
{
  if (not utils::MasterSlave::isSlave()) {

    _iterationsWriter = std::make_shared<io::TXTTableWriter>("precice-" + _localParticipant + "-iterations.log");
    if (not doesFirstStep()) {
      _convergenceWriter = std::make_shared<io::TXTTableWriter>("precice-" + _localParticipant + "-convergence.log");
    }

    // check if coarse model optimization exists
    bool hasCoarseModelOptimization = false;
    for (ConvergenceMeasure &convMeasure : _convergenceMeasures)
      if (convMeasure.level > 0)
        hasCoarseModelOptimization = true;

    _iterationsWriter->addData("Timesteps", io::TXTTableWriter::INT);
    _iterationsWriter->addData("Total_Iterations", io::TXTTableWriter::INT);
    _iterationsWriter->addData("Iterations", io::TXTTableWriter::INT);
    if (hasCoarseModelOptimization) {
      _iterationsWriter->addData("Total_Iterations_Surrogate_Model", io::TXTTableWriter::INT);
      _iterationsWriter->addData("Iterations_Surrogate_Model", io::TXTTableWriter::INT);
    }
    _iterationsWriter->addData("Convergence", io::TXTTableWriter::INT);

    if (not doesFirstStep()) {
      _convergenceWriter->addData("Timestep", io::TXTTableWriter::INT);
      _convergenceWriter->addData("Iteration", io::TXTTableWriter::INT);
    }

    int i = -1;
    if (not doesFirstStep()) {
      for (ConvergenceMeasure &convMeasure : _convergenceMeasures) {
        i++;
        // only for fine model optimization, i.e., coupling
        if (convMeasure.level > 0)
          continue;
        std::stringstream sstm, sstm2;
        sstm << "avgConvRate(" << i << ")";
        sstm2 << "resNorm(" << i << ")";
        _iterationsWriter->addData(sstm.str(), io::TXTTableWriter::DOUBLE);
        _convergenceWriter->addData(sstm2.str(), io::TXTTableWriter::DOUBLE);
      }
      _iterationsWriter->addData("deleted_Columns", io::TXTTableWriter::INT);
    }
  }
}

void BaseCouplingScheme::advanceTXTWriters()
{
  if (not utils::MasterSlave::isSlave()) {

    // check if coarse model optimization exists
    bool hasCoarseModelOptimization = false;
    for (ConvergenceMeasure &convMeasure : _convergenceMeasures)
      if (convMeasure.level > 0)
        hasCoarseModelOptimization = true;

    _iterationsWriter->writeData("Timesteps", _timesteps - 1);
    _iterationsWriter->writeData("Total_Iterations", _totalIterations);
    _iterationsWriter->writeData("Iterations", _iterations);
    if (hasCoarseModelOptimization) {
      _iterationsWriter->writeData("Total_Iterations_Surrogate_Model", _totalIterationsCoarseOptimization);
      _iterationsWriter->writeData("Iterations_Surrogate_Model", _iterationsCoarseOptimization);
    }
    int converged = _iterations < _maxIterations ? 1 : 0;
    _iterationsWriter->writeData("Convergence", converged);

    if (not doesFirstStep()) {
      for (size_t i = 0; i < _convergenceMeasures.size(); i++) {

        // only for fine model optimization, i.e., coupling
        if (_convergenceMeasures[i].level > 0)
          continue;

        std::stringstream sstm;
        sstm << "avgConvRate(" << i << ")";
        if (math::equals(_firstResiduumNorm[i], 0.)) {
          _iterationsWriter->writeData(sstm.str(), std::numeric_limits<double>::infinity());
        } else {
          double avgConvRate = _convergenceMeasures[i].measure->getNormResidual() / _firstResiduumNorm[i];
          _iterationsWriter->writeData(sstm.str(), std::pow(avgConvRate, 1. / (double) _iterations));
        }
      }
      _iterationsWriter->writeData("deleted_Columns", _deletedColumnsPPFiltering);
    }
  }
}

void BaseCouplingScheme::updateTimeAndIterations(
    bool convergence,
    bool convergenceCoarseOptimization)
{
  bool manifoldmapping = false;
  if (getAcceleration().get() != nullptr) {
    manifoldmapping = _acceleration->isMultilevelBasedApproach();
  }

  if (not convergence) {

    // The computed timestep part equals the timestep length, since the
    // timestep remainder is zero. Subtract the timestep length do another
    // coupling iteration.
    PRECICE_ASSERT(math::greater(getComputedTimestepPart(), 0.0));
    _time = _time - _computedTimestepPart;

    // in case of multilevel PP: only increment outer iteration count if surrogate model has converged.
    if (convergenceCoarseOptimization) {
      _totalIterations++;
      _iterations++;
    } else {
      // in case of multilevel PP: increment the iteration count of the surrogate model
      _iterationsCoarseOptimization++;
      _totalIterationsCoarseOptimization++;
    }
  } else {

    _totalIterationsCoarseOptimization++;
    if (not manifoldmapping)
      _totalIterations++;

    _iterationsCoarseOptimization = 1;
    _iterations                   = manifoldmapping ? 0 : 1;
  }
}

void BaseCouplingScheme::timestepCompleted()
{
  PRECICE_TRACE(getTimesteps(), getTime());
  PRECICE_INFO("Timestep completed");
  setIsCouplingTimestepComplete(true);
  setTimesteps(getTimesteps() + 1);
  if (isCouplingOngoing()) {
    PRECICE_DEBUG("Setting require create checkpoint");
    requireAction(constants::actionWriteIterationCheckpoint());
  }
}

bool BaseCouplingScheme::maxIterationsReached()
{
  if (not _isCoarseModelOptimizationActive) {
    return _iterations == _maxIterations;
  } else {
    return _iterationsCoarseOptimization == _maxIterations;
  }
}
}
} // namespace precice, cplscheme
