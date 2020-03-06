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
    std::vector<m2n::PtrM2N>      m2ns,
    constants::TimesteppingMethod dtMethod,
    int                           maxIterations)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, "neverFirstParticipant",
                         localParticipant, localParticipant, m2n::PtrM2N(), maxIterations, Implicit, dtMethod),
      _communications(m2ns)
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  for (size_t i = 0; i < _communications.size(); ++i) {
    DataMap receiveMap;
    DataMap sendMap;
    _receiveDataVector.push_back(receiveMap);
    _sendDataVector.push_back(sendMap);
  }
}

void MultiCouplingScheme::initializeImplicit()
{
  mergeData();                 // merge send and receive data for all pp calls
  setupConvergenceMeasures();  // needs _couplingData configured
  setupDataMatrices(_allData); // Reserve memory and initialize data with zero
  if (getAcceleration().get() != nullptr) {
    PRECICE_CHECK(getAcceleration()->getDataIDs().size() >= 3,
                  "For parallel coupling, the number of coupling data vectors has to be at least 3, not: "
                      << getAcceleration()->getDataIDs().size());
    getAcceleration()->initialize(_allData); // Reserve memory, initialize
  }
}

void MultiCouplingScheme::initializeImplementation()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");

  for (DataMap &dataMap : _sendDataVector) {
    if(anyDataRequiresInitialization(dataMap)) {
      hasToSendInitializedData();
    }
  }
  for (DataMap &dataMap : _receiveDataVector) {
    if(anyDataRequiresInitialization(dataMap)) {
      hasToReceiveInitializedData();
    }
  }
}

void MultiCouplingScheme::exchangeInitialData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");

  if (receivesInitializedData()) {
    receiveData();

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
  if (sendsInitializedData()) {
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
}

void MultiCouplingScheme::explicitAdvance()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
// @todo implement
}

std::pair<bool, bool> MultiCouplingScheme::implicitAdvance()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  PRECICE_DEBUG("Computed full length of iteration");

  receiveData();

  auto designSpecifications = getAcceleration()->getDesignSpecification(_allData);
  bool convergence          = measureConvergence(designSpecifications);

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
    PRECICE_ASSERT(not getIsCoarseModelOptimizationActive());
    m2n->send(getIsCoarseModelOptimizationActive()); //need to do this to match with ParallelCplScheme
  }

  if (convergence && (getExtrapolationOrder() > 0)) {
    extrapolateData(_allData); // Also stores data
  } else {                     // Store data for convergence measurement, acceleration, or extrapolation
    for (DataMap::value_type &pair : _allData) {
      if (pair.second->oldValues.size() > 0) {
        pair.second->oldValues.col(0) = *pair.second->values;
      }
    }
  }

  sendData();

  // TODO: Returning a hard-coded "true" is wrong or at least has a smell! We do not make use of the default value in BaseCouplingScheme::updateTimeAndIterations(bool convergence, bool convergenceCoarseOptimization = true), but provide it explicitly!
  return std::pair<bool, bool>(convergence, true);
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
                            << "added twice for sending.");
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
                            << "added twice for receiving.");
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
  setHasDataBeenExchanged(true);
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

void MultiCouplingScheme::assignDataToConvergenceMeasure(ConvergenceMeasure *convergenceMeasure, int dataID)
{
  convergenceMeasure->couplingData = getData(dataID);
  PRECICE_ASSERT(convergenceMeasure->couplingData != nullptr);
}

} // namespace cplscheme
} // namespace precice
