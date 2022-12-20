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

void MultiCouplingScheme::clearTimeStepSendStorage()
{
  for (auto &sendExchange : _sendDataVector | boost::adaptors::map_values) {
    for (const auto &data : sendExchange | boost::adaptors::map_values) {
      data->clearTimeStepsStorage();
    }
  }
}

void MultiCouplingScheme::storeTimeStepSendData(double relativeDt)
{
  for (auto &sendExchange : _sendDataVector | boost::adaptors::map_values) {
    for (const auto &data : sendExchange | boost::adaptors::map_values) {
      data->storeValuesAtTime(relativeDt, data->values());
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
      }
      checkDataHasBeenReceived();
    } else {
      for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
        initializeZeroReceiveData(receiveExchange);
      }
    }
  }
  PRECICE_DEBUG("Initial data is exchanged in MultiCouplingScheme");
}

void MultiCouplingScheme::storeReceiveData(double relativeDt)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
    for (auto &receiveData : receiveExchange | boost::adaptors::map_values) {
      bool mustOverride = true;
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

typedef std::map<int, PtrCouplingData> DataMap;

const DataMap MultiCouplingScheme::getAllData()
{
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

void MultiCouplingScheme::performReceiveOfFirstAdvance()
{
  return; // no action needed.
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
    // @todo similar code breaks in SerialCouplingScheme.cpp for CplSchemeTests/SerialImplicitCouplingSchemeTests/ testConfiguredAbsConvergenceMeasureSynchronized. Why? @fsimonis
    // for (const auto &data : getAllData() | boost::adaptors::map_values) {
    //   data->moveTimeStepsStorage();
    // }
    for (const DataMap::value_type &data : getAllData()) {
      data.second->moveTimeStepsStorage();
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
