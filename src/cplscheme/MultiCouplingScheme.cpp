#include "MultiCouplingScheme.hpp"
#include <algorithm>
#include <map>
#include <memory>
#include <ostream>
#include <stddef.h>
#include <type_traits>
#include <utility>
#include "acceleration/Acceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"

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
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, localParticipant, maxIterations, Implicit, dtMethod),
      _m2ns(m2ns)
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  setDoesFirstStep(false); // MultiCouplingScheme never does the first step, because it is never the first participant
  for (size_t i = 0; i < _m2ns.size(); ++i) {
    DataMap receiveMap;
    DataMap sendMap;
    _receiveDataVector.push_back(receiveMap);
    _sendDataVector.push_back(sendMap);
  }
}

std::vector<std::string> MultiCouplingScheme::getCouplingPartners() const
{
  std::vector<std::string> partnerNames;
  // Add non-local participant
  // @todo has to be implemented!
  PRECICE_ASSERT(false);
  return partnerNames;
}

void MultiCouplingScheme::initializeImplementation()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");

  for (DataMap &sendData : _sendDataVector) {
    determineInitialSend(sendData);
  }
  for (DataMap &receiveData : _receiveDataVector) {
    determineInitialReceive(receiveData);
  }
}

void MultiCouplingScheme::exchangeInitialData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");

  if (receivesInitializedData()) {
    for (size_t i = 0; i < _m2ns.size(); i++) {
      receiveData(_m2ns[i], _receiveDataVector[i]);
    }
    checkDataHasBeenReceived();
    // second participant has to save values for extrapolation
    for (DataMap &receiveData : _receiveDataVector) {
      updateOldValues(receiveData);
    }
  }
  if (sendsInitializedData()) {
    for (DataMap &sendData : _sendDataVector) {
      updateOldValues(sendData);
    }
    for (size_t i = 0; i < _m2ns.size(); i++) {
      sendData(_m2ns[i], _sendDataVector[i]);
    }
  }
}

bool MultiCouplingScheme::exchangeDataAndAccelerate()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // @todo implement MultiCouplingScheme for explicit coupling

  PRECICE_DEBUG("Computed full length of iteration");

  for (size_t i = 0; i < _m2ns.size(); i++) {
    receiveData(_m2ns[i], _receiveDataVector[i]);
  }
  checkDataHasBeenReceived();

  PRECICE_DEBUG("Perform acceleration (only second participant)...");
  bool convergence = accelerate();

  for (m2n::PtrM2N m2n : _m2ns) {
    sendConvergence(m2n, convergence);
  }

  for (size_t i = 0; i < _m2ns.size(); i++) {
    sendData(_m2ns[i], _sendDataVector[i]);
  }

  return convergence;
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
    PtrCouplingData     ptrCplData(new CouplingData(data, mesh, initialize));
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
    PtrCouplingData     ptrCplData(new CouplingData(data, mesh, initialize));
    DataMap::value_type pair = std::make_pair(id, ptrCplData);
    _receiveDataVector[index].insert(pair);
  } else {
    PRECICE_ERROR("Data \"" << data->getName()
                            << "\" of mesh \"" << mesh->getName() << "\" cannot be "
                            << "added twice for receiving.");
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

void MultiCouplingScheme::assignDataToConvergenceMeasure(ConvergenceMeasureContext *convergenceMeasure, int dataID)
{
  convergenceMeasure->couplingData = getData(dataID);
  PRECICE_ASSERT(convergenceMeasure->couplingData != nullptr);
}

} // namespace cplscheme
} // namespace precice
