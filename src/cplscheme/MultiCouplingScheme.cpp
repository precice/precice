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
  return _allData;
}

void MultiCouplingScheme::overwriteSendValuesAtWindowEnd()
{
  DataMap uniqueSendData;
  for (auto &sendExchange : _sendDataVector | boost::adaptors::map_values) {
    for (auto &pair : sendExchange) {
      uniqueSendData[pair.first] = pair.second;
    }
  }

  for (const auto &data : uniqueSendData | boost::adaptors::map_values) {
    data->overwriteValuesAtWindowEnd(data->values());
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
    data->storeValuesAtTime(time::Storage::WINDOW_START, data->values());
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

void MultiCouplingScheme::storeReceiveData(double relativeDt)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);

  DataMap uniqueReceiveData;
  for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
    for (auto &pair : receiveExchange) {
      uniqueReceiveData[pair.first] = pair.second;
    }
  }

  for (auto &receiveData : uniqueReceiveData | boost::adaptors::map_values) {
    bool mustOverride = true;
    receiveData->storeValuesAtTime(relativeDt, receiveData->values(), mustOverride);
  }
}

void MultiCouplingScheme::loadReceiveDataFromStorage(double relativeDt)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);

  DataMap uniqueReceiveData;
  for (auto &receiveExchange : _receiveDataVector | boost::adaptors::map_values) {
    for (auto &pair : receiveExchange) {
      uniqueReceiveData[pair.first] = pair.second;
    }
  }

  for (auto &receiveData : uniqueReceiveData | boost::adaptors::map_values) {
    receiveData->values() = receiveData->getValuesAtTime(relativeDt);
  }
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
