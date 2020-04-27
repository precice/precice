#include "MultiCouplingScheme.hpp"
#include "acceleration/Acceleration.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "math/math.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {

MultiCouplingScheme::MultiCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    int                           validDigits,
    const std::string &           localParticipant,
    std::vector<m2n::PtrM2N>      m2n,
    constants::TimesteppingMethod dtMethod,
    int                           maxIterations)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, "neverFirstParticipant",
                         localParticipant, localParticipant, m2n::PtrM2N(), maxIterations, dtMethod),
      _communications(m2n)
{
  for (size_t i = 0; i < _communications.size(); ++i) {
    DataMap receiveMap;
    DataMap sendMap;
    _receiveDataVector.push_back(receiveMap);
    _sendDataVector.push_back(sendMap);
  }
}

void MultiCouplingScheme::initialize(
    double startTime,
    int    startTimeWindow)
{
  PRECICE_TRACE(startTime, startTimeWindow);
  PRECICE_ASSERT(not isInitialized());
  PRECICE_ASSERT(math::greaterEquals(startTime, 0.0), startTime);
  PRECICE_ASSERT(startTimeWindow >= 0, startTimeWindow);
  setTime(startTime);
  setTimeWindows(startTimeWindow);

  mergeData();                 // merge send and receive data for all pp calls
  setupConvergenceMeasures();  // needs _couplingData configured
  setupDataMatrices(_allData); // Reserve memory and initialize data with zero
  if (getAcceleration().get() != nullptr) {
    PRECICE_CHECK(getAcceleration()->getDataIDs().size() >= 3,
                  "For parallel coupling, the number of coupling data vectors has to be at least 3, not: "
                      << getAcceleration()->getDataIDs().size());
    getAcceleration()->initialize(_allData); // Reserve memory, initialize
  }

  requireAction(constants::actionWriteIterationCheckpoint());
  initializeTXTWriters();

  for (DataMap &dataMap : _sendDataVector) {
    for (DataMap::value_type &pair : dataMap) {
      if (pair.second->initialize) {
        setHasToSendInitData(true);
        break;
      }
    }
  }
  for (DataMap &dataMap : _receiveDataVector) {
    for (DataMap::value_type &pair : dataMap) {
      if (pair.second->initialize) {
        setHasToReceiveInitData(true);
        break;
      }
    }
  }

  if (hasToSendInitData()) {
    requireAction(constants::actionWriteInitialData());
  }

  setIsInitialized(true);
}

void MultiCouplingScheme::initializeData()
{
  PRECICE_TRACE();
  PRECICE_CHECK(isInitialized(), "initializeData() can be called after initialize() only!");

  if (not hasToSendInitData() && not hasToReceiveInitData()) {
    PRECICE_INFO("initializeData is skipped since no data has to be initialized");
    return;
  }

  PRECICE_CHECK(not(hasToSendInitData() && isActionRequired(constants::actionWriteInitialData())),
                "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);

  if (hasToReceiveInitData()) {
    receiveData();
    setHasDataBeenExchanged(true);

    // second participant has to save values for extrapolation
    if (getExtrapolationOrder() > 0) {
      for (DataMap &dataMap : _receiveDataVector) {
        for (DataMap::value_type &pair : dataMap) {
          pair.second->oldValues.col(0) = *pair.second->values;
          // For extrapolation, treat the initial value as old time window value
          utils::shiftSetFirst(pair.second->oldValues, *pair.second->values);
        }
      }
    }
  }
  if (hasToSendInitData()) {
    if (getExtrapolationOrder() > 0) {
      for (DataMap &dataMap : _sendDataVector) {
        for (DataMap::value_type &pair : dataMap) {
          pair.second->oldValues.col(0) = *pair.second->values;
          // For extrapolation, treat the initial value as old time window value
          utils::shiftSetFirst(pair.second->oldValues, *pair.second->values);
        }
      }
    }
    sendData();
  }

  // in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}

void MultiCouplingScheme::advance()
{
  PRECICE_TRACE(getTimeWindows(), getTime());
  checkCompletenessRequiredActions();

  PRECICE_CHECK(!hasToReceiveInitData() && !hasToSendInitData(),
                "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsTimeWindowComplete(false);
  bool convergence = false;
  if (math::equals(getThisTimeWindowRemainder(), 0.0, _eps)) {
    PRECICE_DEBUG("Computed full length of iteration");

    receiveData();

    convergence = measureConvergence();

    // Stop, when maximal iteration count (given in config) is reached
    if (maxIterationsReached()) {
      convergence = true;
    }
    if (convergence) {
      if (getAcceleration().get() != nullptr) {
        getAcceleration()->iterationsConverged(_allData);
      }
      newConvergenceMeasurements();
      timeWindowCompleted();
    } else if (getAcceleration().get() != nullptr) {
      getAcceleration()->performAcceleration(_allData);
    }

    for (m2n::PtrM2N m2n : _communications) {
      m2n->send(convergence);
    }

    if (convergence && (getExtrapolationOrder() > 0)) {
      extrapolateData(_allData); // Also stores data
    } else {                     // Store data for conv. measurement, acceleration, or extrapolation
      for (DataMap::value_type &pair : _allData) {
        if (pair.second->oldValues.size() > 0) {
          pair.second->oldValues.col(0) = *pair.second->values;
        }
      }
    }
    sendData();

    if (not convergence) {
      PRECICE_DEBUG("No convergence achieved");
      requireAction(constants::actionReadIterationCheckpoint());
    } else {
      PRECICE_DEBUG("Convergence achieved");
      advanceTXTWriters();
    }
    updateTimeAndIterations(convergence);
    setHasDataBeenExchanged(true);
    setComputedTimeWindowPart(0.0);
  } // subcycling complete
}

void MultiCouplingScheme::mergeData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_allData.empty(), "This function should only be called once.");
  PRECICE_ASSERT(_sendDataVector.size() == _receiveDataVector.size());
  for (size_t i = 0; i < _sendDataVector.size(); i++) {
    _allData.insert(_sendDataVector[i].begin(), _sendDataVector[i].end());
    _allData.insert(_receiveDataVector[i].begin(), _receiveDataVector[i].end());
  }
}

void MultiCouplingScheme::addDataToSend(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize,
    int           index)
{
  int id = data->getID();
  if (!utils::contained(id, _sendDataVector[index])) {
    PtrCouplingData     ptrCplData(new CouplingData(&(data->values()), mesh, initialize, data->getDimensions()));
    DataMap::value_type pair = std::make_pair(id, ptrCplData);
    _sendDataVector[index].insert(pair);
  } else {
    PRECICE_ERROR("Data \"" << data->getName()
                            << "\" of mesh \"" << mesh->getName() << "\" cannot be "
                            << "added twice for sending!");
  }
}

void MultiCouplingScheme::addDataToReceive(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize,
    int           index)
{
  int id = data->getID();
  if (!utils::contained(id, _receiveDataVector[index])) {
    PtrCouplingData     ptrCplData(new CouplingData(&(data->values()), mesh, initialize, data->getDimensions()));
    DataMap::value_type pair = std::make_pair(id, ptrCplData);
    _receiveDataVector[index].insert(pair);
  } else {
    PRECICE_ERROR("Data \"" << data->getName()
                            << "\" of mesh \"" << mesh->getName() << "\" cannot be "
                            << "added twice for receiving!");
  }
}

void MultiCouplingScheme::sendData()
{
  PRECICE_TRACE();

  for (size_t i = 0; i < _communications.size(); i++) {
    PRECICE_ASSERT(_communications[i].get() != nullptr);
    PRECICE_ASSERT(_communications[i]->isConnected());

    for (DataMap::value_type &pair : _sendDataVector[i]) {
      int size = pair.second->values->size();
      if (size > 0) {
        _communications[i]->send(pair.second->values->data(), size, pair.second->mesh->getID(), pair.second->dimension);
      }
    }
  }
}

void MultiCouplingScheme::receiveData()
{
  PRECICE_TRACE();

  for (size_t i = 0; i < _communications.size(); i++) {
    PRECICE_ASSERT(_communications[i].get() != nullptr);
    PRECICE_ASSERT(_communications[i]->isConnected());

    for (DataMap::value_type &pair : _receiveDataVector[i]) {
      int size = pair.second->values->size();
      if (size > 0) {
        _communications[i]->receive(pair.second->values->data(), size, pair.second->mesh->getID(), pair.second->dimension);
      }
    }
  }
}

void MultiCouplingScheme::setupConvergenceMeasures()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not doesFirstStep());
  PRECICE_ASSERT(not _convergenceMeasures.empty(),
                 "At least one convergence measure has to be defined for an implicit coupling scheme!");
  for (ConvergenceMeasure &convMeasure : _convergenceMeasures) {
    int dataID               = convMeasure.data->getID();
    convMeasure.couplingData = getData(dataID);
    PRECICE_ASSERT(convMeasure.couplingData != nullptr);
  }
}

CouplingData *MultiCouplingScheme::getData(
    int dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _allData.find(dataID);
  if (iter != _allData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

} // namespace cplscheme
} // namespace precice
