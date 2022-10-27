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
#include "com/Communication.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "utils/IntraComm.hpp"

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

CouplingScheme::ChangedMeshes MultiCouplingScheme::firstSynchronization(const CouplingScheme::ChangedMeshes &changes)
{
  PRECICE_DEBUG("Exchanging mesh changes...");
  if (_isController) {
    sendLocalChanges(changes);
    return receiveRemoteChanges();
  } else {
    auto remote = receiveRemoteChanges();
    sendLocalChanges(changes);
    return remote;
  }
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

CouplingScheme::ChangedMeshes MultiCouplingScheme::secondSynchronization()
{
  return {};
}

void MultiCouplingScheme::exchangeSecondData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // @todo implement MultiCouplingScheme for explicit coupling

  if (_isController) {
    doImplicitStep();
    for (const auto &m2nPair : _m2ns) {
      sendConvergence(m2nPair.second);
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

CouplingScheme::ChangedMeshes MultiCouplingScheme::receiveRemoteChanges()
{
  CouplingScheme::ChangedMeshes changes;
  /// Receive to primary rank
  if (_isController) {
    if (!utils::IntraComm::isSecondary()) {
      // Receive changes from all non-controllers
      std::set<MeshID> allChanges;
      for (auto &m2n : _m2ns) {
        auto changes = m2n.second->getPrimaryRankCommunication()->receiveRange(0, com::AsVectorTag<int>{});
        allChanges.insert(changes.begin(), changes.end());
      }
      // Broadcast to other secondaries of the controller
      changes.insert(changes.end(), allChanges.begin(), allChanges.end());
    }
  } else {
    if (!utils::IntraComm::isSecondary()) {
      // Receive changes from the controller
      CouplingScheme::ChangedMeshes changes = _m2ns[_controller]->getPrimaryRankCommunication()->receiveRange(0, com::AsVectorTag<int>{});
    }
  }

  /// Broadcast changes to secondaries rank
  if (utils::IntraComm::isPrimary()) {
    // Broadcast to other secondaries of the controller
    utils::IntraComm::getCommunication()->broadcast(changes);
  }
  if (utils::IntraComm::isSecondary()) {
    // Broadcast to other secondaries of the controller
    utils::IntraComm::getCommunication()->broadcast(changes, 0);
  }
  return changes;
}

void MultiCouplingScheme::sendLocalChanges(const CouplingScheme::ChangedMeshes &changes)
{
  // Gather changes on primary of this participant
  if (utils::IntraComm::isSecondary()) {
    utils::IntraComm::getCommunication()->sendRange(changes, 0);
    return;
  }

  std::set<MeshID> allChanges{changes.begin(), changes.end()};
  if (utils::IntraComm::isParallel()) {
    for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
      auto next = utils::IntraComm::getCommunication()->receiveRange(rank, com::AsVectorTag<int>{});
      allChanges.insert(next.begin(), next.end());
    }
  }
  CouplingScheme::ChangedMeshes changesToSend(allChanges.begin(), allChanges.end());

  // Send changes
  if (_isController) {
    for (auto &m2n : _m2ns) {
      m2n.second->getPrimaryRankCommunication()->sendRange(changesToSend, 0);
    }
  } else {
    _m2ns[_controller]->getPrimaryRankCommunication()->sendRange(changesToSend, 0);
  }
}

} // namespace precice::cplscheme
