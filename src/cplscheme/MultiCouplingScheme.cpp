#include "MultiCouplingScheme.hpp"
#include <algorithm>
#include <boost/range/adaptor/map.hpp>
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
  for (auto &sendExchange : _sendDataVector | boost::adaptors::map_values) {
    determineInitialSend(sendExchange);
  }
  for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
    determineInitialReceive(receiveExchange);
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

bool MultiCouplingScheme::hasAnySendData()
{
  return std::any_of(_sendDataVector.cbegin(), _sendDataVector.cend(), [](const auto &sendExchange) { return not sendExchange.second.empty(); });
}

const DataMap MultiCouplingScheme::getAccelerationData()
{
  // MultiCouplingScheme applies acceleration to all CouplingData
  return _allData;
}

void MultiCouplingScheme::storeSendValuesAtTime(double relativeDt)
{
  DataMap uniqueSendData;
  for (auto &sendExchange : _sendDataVector | boost::adaptors::map_values) {
    for (auto &pair : sendExchange) {
      uniqueSendData[pair.first] = pair.second;
    }
  }

  for (auto &data : uniqueSendData | boost::adaptors::map_values) {
    data->storeValuesAtTime(relativeDt, data->values());
  }
}

void MultiCouplingScheme::initializeSendDataStorage()
{
  DataMap uniqueSendData;
  for (auto &sendExchange : _sendDataVector | boost::adaptors::map_values) {
    for (auto &pair : sendExchange) {
      uniqueSendData[pair.first] = pair.second;
    }
  }

  for (const auto &data : uniqueSendData | boost::adaptors::map_values) {
    // initialize as constant
    data->storeValuesAtTime(time::Storage::WINDOW_START, data->values());
    data->storeValuesAtTime(time::Storage::WINDOW_END, data->values());
  }
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
    } else {
      DataMap uniqueReceiveData;
      for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
        for (auto &pair : receiveExchange) {
          uniqueReceiveData[pair.first] = pair.second;
        }
      }
      initializeZeroReceiveData(uniqueReceiveData);
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
    } else {
      DataMap uniqueReceiveData;
      for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
        for (auto &pair : receiveExchange) {
          uniqueReceiveData[pair.first] = pair.second;
        }
      }
      initializeZeroReceiveData(uniqueReceiveData);
    }
  }
  PRECICE_DEBUG("Initial data is exchanged in MultiCouplingScheme");
}

bool MultiCouplingScheme::hasReceiveData(std::string dataName)
{
  for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
    for (auto &receiveData : receiveExchange | boost::adaptors::map_values) {
      if (receiveData->getDataName() == dataName) {
        return true;
      }
    }
  }
  return false;
}

void MultiCouplingScheme::loadReceiveDataFromStorage(std::string dataName, double relativeDt)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  // @todo work with _allData and move into BaseCouplingScheme
  for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
    for (auto &receiveData : receiveExchange | boost::adaptors::map_values) {
      if (receiveData->getDataName() == dataName) {
        receiveData->values() = receiveData->getValuesAtTime(relativeDt);
        return;
      }
    }
  }
  PRECICE_ASSERT(false, "Data with name not found", dataName);
}

// @todo may be moved into BaseCouplingScheme, but should be done consistently with MultiCouplingScheme::loadReceiveDataFromStorage
void MultiCouplingScheme::clearAllDataStorage()
{
  for (auto &data : _allData | boost::adaptors::map_values) {
    data->clearTimeStepsStorage();
  }
}

std::vector<double> MultiCouplingScheme::getReceiveTimes(std::string dataName)
{
  auto times = std::vector<double>();
  for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
    for (auto &data : receiveExchange | boost::adaptors::map_values) {
      if (data->getDataName() == dataName) {
        auto timesVec = data->getStoredTimesAscending();
        PRECICE_ASSERT(timesVec.size() > 0, timesVec.size());
        for (int i = 0; i < timesVec.size(); i++) {
          times.push_back(timesVec(i));
        }
        PRECICE_DEBUG("Receive times for {} are {}", dataName, times);
        return times;
      }
    }
  }
  PRECICE_DEBUG("No data with dataName {} found in receive data. Returning empty.", dataName);
  PRECICE_ASSERT(false);
  return times;
}

void MultiCouplingScheme::exchangeFirstData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // @todo implement MultiCouplingScheme for explicit coupling

  PRECICE_DEBUG("Computed full length of iteration");

  if (_isController) {
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();

  } else {
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  }
}

void MultiCouplingScheme::exchangeSecondData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // @todo implement MultiCouplingScheme for explicit coupling

  if (_isController) {
    PRECICE_DEBUG("Perform acceleration (only controller)...");
    doImplicitStep();
    for (const auto &m2n : _m2ns | boost::adaptors::map_values) {
      sendConvergence(m2n);
    }
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  } else {
    receiveConvergence(_m2ns[_controller]);
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
  }

  if (hasConverged()) {
    for (const auto &data : _allData | boost::adaptors::map_values) {
      data->moveTimeStepsStorage();
    }
  }
  if (isImplicitCouplingScheme()) {
    storeIteration();
  }
}

void MultiCouplingScheme::addDataToSend(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 requiresInitialization,
    const std::string &  to)
{
  PtrCouplingData ptrCplData = addCouplingData(data, std::move(mesh), requiresInitialization);
  PRECICE_DEBUG("Configuring send data to {}", to);
  _sendDataVector[to].emplace(data->getID(), ptrCplData);
}

void MultiCouplingScheme::addDataToReceive(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 requiresInitialization,
    const std::string &  from)
{
  PtrCouplingData ptrCplData = addCouplingData(data, std::move(mesh), requiresInitialization);
  PRECICE_DEBUG("Configuring receive data from {}", from);
  _receiveDataVector[from].emplace(data->getID(), ptrCplData);
}

} // namespace precice::cplscheme
