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
    double                             maxTime,
    int                                maxTimeWindows,
    double                             timeWindowSize,
    int                                validDigits,
    const std::string &                localParticipant,
    std::map<std::string, m2n::PtrM2N> m2ns,
    constants::TimesteppingMethod      dtMethod,
    std::string                        controller,
    int                                maxIterations)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, localParticipant, maxIterations, Implicit, dtMethod),
      _m2ns(std::move(m2ns)), _isController(localParticipant == controller)
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");

  PRECICE_DEBUG("CONTROLLER IS {}", controller);

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
  for (auto &sendExchange : _sendDataVector) {
    determineInitialSend(sendExchange.second);
  }
  for (auto &receiveExchange : _receiveDataVector) {
    determineInitialReceive(receiveExchange.second);
  }
  PRECICE_DEBUG("MultiCouplingScheme is initialized.");
}

void MultiCouplingScheme::exchangeInitialData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");

  if (_isController) {
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
      }

      checkDataHasBeenReceived();
      // second participant has to save values for extrapolation
      for (auto &receiveExchange : _receiveDataVector) {
        updateOldValues(receiveExchange.second);
      }
    }
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        updateOldValues(sendExchange.second);
      }
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
  } else {
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        updateOldValues(sendExchange.second);
      }
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
      }
      checkDataHasBeenReceived();
      // second participant has to save values for extrapolation
      for (auto &receiveExchange : _receiveDataVector) {
        updateOldValues(receiveExchange.second);
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

  if (_isController) {
    for (auto &receiveExchange : _receiveDataVector) {
      PRECICE_DEBUG("REQUESTED RECEIVE IS: {}", receiveExchange.first);
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();

    PRECICE_DEBUG("Perform acceleration (only second participant)...");
    convergence = accelerate();

    for (const auto &m2nPair : _m2ns) {
      sendConvergence(m2nPair.second, convergence);
    }
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  } else {
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }

    bool localConvergence = true;

    for (const auto &m2nPair : _m2ns) {
      localConvergence = receiveConvergence(m2nPair.second);
      convergence      = (convergence && localConvergence);
    }

    for (auto &receiveExchange : _receiveDataVector) {
      PRECICE_DEBUG("REQUESTED RECEIVE IS: {}", receiveExchange.first);
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
  }

  return convergence;
}

void MultiCouplingScheme::mergeData()
{
  // PRECICE_TRACE();
  // PRECICE_ASSERT(_allData.empty(), "This function should only be called once.");
  // PRECICE_ASSERT(_sendDataVector.size() == _receiveDataVector.size());
}

void MultiCouplingScheme::addDataToSend(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize,
    std::string   to)
{
  int id = data->getID();
  PRECICE_DEBUG("Configuring send data to {}", to);
  PtrCouplingData     ptrCplData(new CouplingData(data, mesh, initialize));
  DataMap::value_type dataPair = std::make_pair(id, ptrCplData);
  _sendDataVector[to].insert(dataPair);
  _allData.insert(dataPair);
}

void MultiCouplingScheme::addDataToReceive(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize,
    std::string   from)
{
  int id = data->getID();
  PRECICE_DEBUG("Configuring receive data from {}", from);
  PtrCouplingData     ptrCplData(new CouplingData(data, mesh, initialize));
  DataMap::value_type dataPair = std::make_pair(id, ptrCplData);
  _receiveDataVector[from].insert(dataPair);
  _allData.insert(dataPair);
}

CouplingData *MultiCouplingScheme::getData(
    int dataID)
{
  PRECICE_TRACE(dataID);
  PRECICE_DEBUG("Desired data ID is {}", dataID);
  DataMap::iterator iter = _allData.find(dataID);
  if (iter != _allData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

void MultiCouplingScheme::assignDataToConvergenceMeasure(ConvergenceMeasureContext *convergenceMeasure, int dataID)
{
  //if(_isController){
  convergenceMeasure->couplingData = getData(dataID);
  PRECICE_ASSERT(convergenceMeasure->couplingData != nullptr);
  //}
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
