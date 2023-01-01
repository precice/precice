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
  return getAllData();
}

void MultiCouplingScheme::storeSendValuesAtTime(double relativeDt)
{
  for (auto &sendExchange : _sendDataVector | boost::adaptors::map_values) {
    for (auto &data : sendExchange | boost::adaptors::map_values) {
      data->storeValuesAtTime(relativeDt, data->values());
    }
  }
}

void MultiCouplingScheme::initializeSendDataStorage()
{
  for (auto &sendExchange : _sendDataVector | boost::adaptors::map_values) {
    for (const auto &data : sendExchange | boost::adaptors::map_values) {
      data->storeValuesAtTime(time::Storage::WINDOW_START, data->values());
    }
  }
}

void MultiCouplingScheme::exchangeInitialData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  bool initialCommunication = true;

  if (_isController) {
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second, initialCommunication);
        storeReceiveData(time::Storage::WINDOW_END);
      }
      checkDataHasBeenReceived();
    } else {
      for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
        initializeZeroReceiveData(receiveExchange);
      }
    }
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second, initialCommunication);
      }
    }
  } else {
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second, initialCommunication);
      }
    }
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second, initialCommunication);
        storeReceiveData(time::Storage::WINDOW_END);
      }
      checkDataHasBeenReceived();
    } else {
      for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
        initializeZeroReceiveData(receiveExchange);
      }
    }
  }

  for (auto &sendExchange : _sendDataVector) {
    for (const auto &data : sendExchange.second | boost::adaptors::map_values) {
      data->clearTimeStepsStorage();
    }
  }
  PRECICE_DEBUG("Initial data is exchanged in MultiCouplingScheme");
}

void MultiCouplingScheme::storeReceiveData(double relativeDt, bool mustOverride)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
    for (auto &receiveData : receiveExchange | boost::adaptors::map_values) {
      receiveData->storeValuesAtTime(relativeDt, receiveData->values(), mustOverride);
    }
  }
}

void MultiCouplingScheme::loadReceiveDataFromStorage(double relativeDt)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
    for (auto &receiveData : receiveExchange | boost::adaptors::map_values) {
      receiveData->values() = receiveData->getValuesAtTime(relativeDt);
    }
  }
}

std::vector<double> MultiCouplingScheme::getReceiveTimes()
{
  //@todo stub implementation. Should walk over all receive data, get times and ensure that all times vectors actually hold the same times (since otherwise we would have to get times individually per data)
  //@todo subcycling is not supported for MultiCouplingScheme, because this needs a complicated interplay of picking the right data in time and mapping this data. This is hard to realize with the current implementation.
  auto times = std::vector<double>({time::Storage::WINDOW_END});
  return times;
}

const DataMap MultiCouplingScheme::getAllData()
{
  // MultiCouplingScheme has to collect all send data and receive data from _sendDataVector and _receiveDataVector
  DataMap allData;
  // @todo use C++17 std::map::merge
  for (auto &sendData : _sendDataVector | boost::adaptors::map_values) {
    allData.insert(sendData.begin(), sendData.end());
  }
  for (auto &receiveData : _receiveDataVector | boost::adaptors::map_values) {
    allData.insert(receiveData.begin(), receiveData.end());
  }
  return allData;
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
  } else {
    receiveConvergence(_m2ns[_controller]);
  }

  if (hasConverged()) {
    // @todo similar code breaks in SerialCouplingScheme.cpp for CplSchemeTests/SerialImplicitCouplingSchemeTests/ testConfiguredAbsConvergenceMeasureSynchronized. Why? @fsimonis
    // for (const auto &data : getAllData() | boost::adaptors::map_values) {
    //   data->moveTimeStepsStorage();
    // }
    for (const DataMap::value_type &data : getAllData()) {
      data.second->moveTimeStepsStorage();
    }
  }
  if (isImplicitCouplingScheme()) {
    storeIteration();
  }

  if (_isController) {
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  } else {
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
  }

  for (auto &sendExchange : _sendDataVector | boost::adaptors::map_values) {
    for (auto &sendData : sendExchange | boost::adaptors::map_values) {
      sendData->clearTimeStepsStorage();
    }
  }
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
