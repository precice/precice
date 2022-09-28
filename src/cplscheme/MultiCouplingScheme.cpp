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

namespace precice::cplscheme {

MultiCouplingScheme::MultiCouplingScheme(
    double                             maxTime,
    int                                maxTimeWindows,
    double                             timeWindowSize,
    int                                validDigits,
    const std::string &                localParticipant,
    std::map<std::string, m2n::PtrM2N> m2ns,
    constants::TimesteppingMethod      dtMethod,
    const std::string &                controller,
    int                                maxIterations,
    int                                extrapolationOrder)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, localParticipant, maxIterations, Implicit, dtMethod, extrapolationOrder),
      _m2ns(std::move(m2ns)), _controller(controller), _isController(controller == localParticipant)
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // Controller participant never does the first step, because it is never the first participant
  setDoesFirstStep(!_isController);
  PRECICE_DEBUG("MultiCoupling scheme is created for {}.", localParticipant);
}

void MultiCouplingScheme::determineInitialDataExchange()
{
  for (auto &sendExchange : _sendDataVector) {
    determineInitialSend(sendExchange.second);
  }
  for (auto &receiveExchange : _receiveDataVector) {
    determineInitialReceive(receiveExchange.second);
  }
}

std::vector<std::string> MultiCouplingScheme::getCouplingPartners() const
{
  std::vector<std::string> partnerNames;

  for (const auto &m2nPair : _m2ns) {
    partnerNames.emplace_back(m2nPair.first);
  }

  return partnerNames;
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
    }
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
  } else {
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
      }
      checkDataHasBeenReceived();
    }
  }
  PRECICE_DEBUG("Initial data is exchanged in MultiCouplingScheme");
}

void MultiCouplingScheme::storeTimeStepSendData(double relativeDt)
{
  PRECICE_ASSERT(relativeDt > 0);
  PRECICE_ASSERT(relativeDt <= 1.0);
  for (auto &sendDataVector : _sendDataVector) {
    for (auto &aSendData : sendDataVector.second) {
      auto theData = aSendData.second->values();
      aSendData.second->storeDataAtTime(theData, relativeDt);
    }
  }
}

void MultiCouplingScheme::storeTimeStepReceiveDataEndOfWindow()
{
  if (hasDataBeenReceived()) {
    for (auto &receiveDataVector : _receiveDataVector) {
      for (auto &aReceiveData : receiveDataVector.second) {
        auto theData = aReceiveData.second->values();
        aReceiveData.second->overrideDataAtEndWindowTime(theData);
      }
    }
  }
}

void MultiCouplingScheme::retreiveTimeStepReceiveData(double relativeDt)
{
  // @todo breaks, if different receiveData live on different time-meshes. This is a realistic use-case for multi coupling! Should use a different signature here to individually retreiveTimeStepData from receiveData. Would also be helpful for mapping.
  PRECICE_ASSERT(relativeDt > 0);
  PRECICE_ASSERT(relativeDt <= 1.0, relativeDt);
  for (auto &receiveDataVector : _receiveDataVector) {
    for (auto &aReceiveData : receiveDataVector.second) {
      aReceiveData.second->values() = aReceiveData.second->getDataAtTime(relativeDt);
      //retreiveTimeStepForData(relativeDt, aReceiveData.second->getDataID());  // @todo: Causes problems. Why? Because DataID is not unique for MultiCouplingScheme's allData.
    }
  }
}

std::vector<double> MultiCouplingScheme::getReceiveTimes()
{
  //@todo stub implementation. Should walk over all receive data, get times and ensure that all times vectors actually hold the same times (since otherwise we would have to get times individually per data)
  //@todo subcycling is not supported for MultiCouplingScheme, because this needs a complicated interplay of picking the right data in time and mapping this data. This is hard to realize with the current implementation.
  auto times = std::vector<double>({1.0});
  return times;
}

typedef std::map<int, PtrCouplingData> DataMap;

const DataMap MultiCouplingScheme::getAllData()
{
  DataMap allData;
  PRECICE_INFO("##### assembling allData() #####");
  // @todo user C++17 std::map::merge
  for (auto &sendData : _sendDataVector) {
    for (auto &aSendData : sendData.second) {
      PRECICE_INFO("DataID: {} {}", aSendData.first, aSendData.second->getDataID());
    }
    allData.insert(sendData.second.begin(), sendData.second.end());
  }
  for (auto &receiveData : _receiveDataVector) {
    for (auto &aReceiveData : receiveData.second) {
      PRECICE_INFO("DataID: {} {}", aReceiveData.first, aReceiveData.second->getDataID());
    }
    allData.insert(receiveData.second.begin(), receiveData.second.end());
  }
  PRECICE_INFO("##### assembling allData() done #####");

  PRECICE_INFO("##### checking allData() #####");
  for (auto &aData : allData) {
    PRECICE_INFO("DataID: {} {}", aData.first, aData.second->getDataID());
  }
  PRECICE_INFO("##### checking allData() done #####");
  return allData;
}

bool MultiCouplingScheme::exchangeDataAndAccelerate()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // @todo implement MultiCouplingScheme for explicit coupling

  PRECICE_DEBUG("Computed full length of iteration");

  bool convergence = true;

  if (_isController) {
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();

    convergence = doImplicitStep();
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

    convergence = receiveConvergence(_m2ns[_controller]);

    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
  }
  return convergence;
}

void MultiCouplingScheme::addDataToSend(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 initialize,
    const std::string &  to)
{
  int id = data->getID();
  PRECICE_DEBUG("Configuring send data to {}", to);
  PtrCouplingData ptrCplData(new CouplingData(data, std::move(mesh), initialize, getExtrapolationOrder()));
  _sendDataVector[to].emplace(id, ptrCplData);
}

void MultiCouplingScheme::addDataToReceive(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 initialize,
    const std::string &  from)
{
  int id = data->getID();
  PRECICE_DEBUG("Configuring receive data from {}", from);
  PtrCouplingData ptrCplData(new CouplingData(data, std::move(mesh), initialize, getExtrapolationOrder()));
  _receiveDataVector[from].emplace(id, ptrCplData);
}

} // namespace precice::cplscheme
