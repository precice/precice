#include "MultiCouplingScheme.hpp"
#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <ostream>
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
    std::map<std::string, m2n::PtrM2N>      m2ns,
    constants::TimesteppingMethod dtMethod,
    int                           maxIterations,
    bool                          isController)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, localParticipant, maxIterations, Implicit, dtMethod),
      _m2ns(std::move(m2ns)), _isController(isController)
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
   
  // Controller participant never does the first step, because it is never the first participant
  setDoesFirstStep(!_isController);

  PRECICE_DEBUG("MultiCoupling scheme is created for {}.", localParticipant);
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

  PRECICE_DEBUG("MultiCouplingScheme is being initialized.");
  for (auto &sendData : _sendDataVector) {
    determineInitialSend(sendData.data);
  }
  for (auto &receiveData : _receiveDataVector) {
    determineInitialReceive(receiveData.data);
  }
  PRECICE_DEBUG("MultiCouplingScheme is initialized.");
}

void MultiCouplingScheme::exchangeInitialData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");

  if(_isController){
    if (receivesInitializedData()) {
      for(auto& receiveExchange : _receiveDataVector){
        receiveData(_m2ns[receiveExchange.with], receiveExchange.data);
      }

      checkDataHasBeenReceived();
      // second participant has to save values for extrapolation
      for (auto &receiveExchange : _receiveDataVector) {
        updateOldValues(receiveExchange.data);
      }
    }
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        updateOldValues(sendExchange.data);
      }
      for (auto& sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.with], sendExchange.data);
      }
    }
  } else {
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        updateOldValues(sendExchange.data);
      }
      for (auto& sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.with], sendExchange.data);
      }
    }
    if (receivesInitializedData()) {
      for(auto& receiveExchange : _receiveDataVector){
        receiveData(_m2ns[receiveExchange.with], receiveExchange.data);
      }
      checkDataHasBeenReceived();
      // second participant has to save values for extrapolation
      for (auto &receiveExchange : _receiveDataVector) {
        updateOldValues(receiveExchange.data);
      }
    }
  }
  PRECICE_DEBUG("Initial data is exchanged in MultiCouplingScheme");
}

bool MultiCouplingScheme::exchangeDataAndAccelerate()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // @todo implement MultiCouplingScheme for explicit coupling

  bool convergence = true;

  PRECICE_DEBUG("Computed full length of iteration");

  if(_isController){
    for(auto& receiveExchange : _receiveDataVector){
      receiveData(_m2ns[receiveExchange.with], receiveExchange.data);
    }
    checkDataHasBeenReceived();

    PRECICE_DEBUG("Perform acceleration (only second participant)...");
    convergence = accelerate();

    for (const auto &m2nPair : _m2ns) {
      sendConvergence(m2nPair.second, convergence);
    }
    for (auto& sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.with], sendExchange.data);
    }
  } else {
    for (auto& sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.with], sendExchange.data);
    }

    bool localConvergence = true;

    for (const auto &m2nPair : _m2ns) {
      localConvergence = receiveConvergence(m2nPair.second);
      convergence = (convergence && localConvergence);
    }
    
    for(auto& receiveExchange : _receiveDataVector){
      receiveData(_m2ns[receiveExchange.with], receiveExchange.data);
    }
    checkDataHasBeenReceived();
  }

  return convergence;
}

void MultiCouplingScheme::mergeData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_allData.empty(), "This function should only be called once.");
  PRECICE_ASSERT(_sendDataVector.size() == _receiveDataVector.size());

  for(size_t i = 0; i < _sendDataVector.size(); ++i){  
    _allData.insert(_sendDataVector[i].data.cbegin(), _sendDataVector[i].data.cend());
    _allData.insert(_receiveDataVector[i].data.cbegin(), _receiveDataVector[i].data.cend());
  }
}

void MultiCouplingScheme::addDataToSend(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize,
    std::string   to)
{
  int id = data->getID();
  if (!utils::contained(id, (*std::find(_sendDataVector.cbegin(), _sendDataVector.cend(), to)).data)) {
    PRECICE_DEBUG("Configuring send data to {}", to);
    PtrCouplingData     ptrCplData(new CouplingData(data, mesh, initialize));
    DataMap dataMap   = {{id, std::move(ptrCplData)}};
    _sendDataVector.emplace_back(std::move(dataMap), to);
  } else {
    PRECICE_ERROR("Data \"{}\" of mesh \"{}\" cannot be added twice for sending.",
                  data->getName(), mesh->getName());
  }
}

void MultiCouplingScheme::addDataToReceive(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize,
    std::string   from)
{
  int id = data->getID(); 
  if (!utils::contained(id, (*std::find(_receiveDataVector.cbegin(), _receiveDataVector.cend(), from)).data)) {
    PRECICE_DEBUG("Configuring receive data from {}", from);
    PtrCouplingData     ptrCplData(new CouplingData(data, mesh, initialize));
    DataMap dataMap = {{id, std::move(ptrCplData)}};
    _receiveDataVector.emplace_back(std::move(dataMap), from);
  } else {
    PRECICE_ERROR("Data \"{}\" of mesh \"{}\" cannot be added twice for receiving.",
                  data->getName(), mesh->getName());
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

bool MultiCouplingScheme::receiveConvergence(m2n::PtrM2N m2n)
{
  PRECICE_ASSERT((!_isController), "For convergence information the receiving participant is always the non-controller one.");
  bool convergence;
  m2n->receive(convergence);
  return convergence;
}

} // namespace cplscheme
} // namespace precice
