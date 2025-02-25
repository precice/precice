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
    const std::string &                localParticipant,
    std::map<std::string, m2n::PtrM2N> m2ns,
    const std::string &                controller,
    int                                minIterations,
    int                                maxIterations)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, localParticipant, minIterations, maxIterations, Implicit, constants::TimesteppingMethod::FIXED_TIME_WINDOW_SIZE),
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

DataMap &MultiCouplingScheme::getAccelerationData()
{
  // MultiCouplingScheme applies acceleration to all CouplingData
  return _allData;
}

void MultiCouplingScheme::initializeReceiveDataStorage()
{
  // @todo check receiveData. Should only contain zero data!
  for (auto &receiveExchange : _receiveDataVector) {
    initializeWithZeroInitialData(receiveExchange.second);
  }
}

void MultiCouplingScheme::exchangeInitialData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  PRECICE_ASSERT(math::equals(getTime(), getWindowStartTime()), getTime(), getWindowStartTime());

  if (_isController) {
    for (auto &[from, data] : _receiveDataVector) {
      if (receivesInitializedDataFrom(from)) {
        receiveData(_m2ns.at(from), data);
      } else {
        initializeWithZeroInitialData(data);
      }
    }
    notifyDataHasBeenReceived();
    for (auto &[to, data] : _sendDataVector) {
      if (sendsInitializedDataTo(to)) {
        sendData(_m2ns.at(to), data);
      }
    }
  } else {
    for (auto &[to, data] : _sendDataVector) {
      if (sendsInitializedDataTo(to)) {
        sendData(_m2ns.at(to), data);
      }
    }
    for (auto &[from, data] : _receiveDataVector) {
      if (receivesInitializedDataFrom(from)) {
        receiveData(_m2ns.at(from), data);
      } else {
        initializeWithZeroInitialData(data);
      }
    }
    notifyDataHasBeenReceived();
  }
  PRECICE_DEBUG("Initial data is exchanged in MultiCouplingScheme");
}

bool MultiCouplingScheme::sendsInitializedDataTo(const std::string &to) const
{
  return _sendInitialTo.count(to) > 0;
}

bool MultiCouplingScheme::receivesInitializedDataFrom(const std::string &from) const
{
  return _receiveInitialFrom.count(from) > 0;
}

void MultiCouplingScheme::exchangeFirstData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  PRECICE_ASSERT(math::equals(getTime(), getWindowEndTime()), getTime(), getWindowEndTime());
  // @todo implement MultiCouplingScheme for explicit coupling

  PRECICE_DEBUG("Computed full length of iteration");

  if (_isController) {
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    notifyDataHasBeenReceived();
  } else {
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  }
}

void MultiCouplingScheme::exchangeSecondData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  PRECICE_ASSERT(math::equals(getTime(), getWindowEndTime()), getTime(), getWindowEndTime());
  // @todo implement MultiCouplingScheme for explicit coupling

  if (not _isController) {
    receiveConvergence(_m2ns[_controller]);
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    notifyDataHasBeenReceived();
  } else {
    doImplicitStep();
    for (const auto &m2n : _m2ns | boost::adaptors::map_values) {
      sendConvergence(m2n);
    }
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  }

  if (hasConverged()) {
    moveToNextWindow();
  }

  storeIteration();
}

void MultiCouplingScheme::addDataToSend(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 requiresInitialization,
    bool                 exchangeSubsteps,
    const std::string &  to)
{
  PtrCouplingData ptrCplData = addCouplingData(data, std::move(mesh), requiresInitialization, exchangeSubsteps, CouplingData::Direction::Send);
  PRECICE_DEBUG("Configuring send data to {}", to);
  _sendDataVector[to].emplace(data->getID(), ptrCplData);
  if (requiresInitialization) {
    _sendInitialTo.emplace(to);
  }
}

void MultiCouplingScheme::addDataToReceive(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 requiresInitialization,
    bool                 exchangeSubsteps,
    const std::string &  from)
{
  PtrCouplingData ptrCplData = addCouplingData(data, std::move(mesh), requiresInitialization, exchangeSubsteps, CouplingData::Direction::Receive);
  PRECICE_DEBUG("Configuring receive data from {}", from);
  _receiveDataVector[from].emplace(data->getID(), ptrCplData);
  if (requiresInitialization) {
    _receiveInitialFrom.emplace(from);
  }
}

} // namespace precice::cplscheme
