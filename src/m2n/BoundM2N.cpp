#include <memory>

#include "com/Communication.hpp"
#include "com/Extra.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/BoundM2N.hpp"
#include "m2n/M2N.hpp"
#include "precice/impl/Types.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"
#include "profiling/Event.hpp"

using precice::profiling::Event;

namespace precice::m2n {

void BoundM2N::prepareEstablishment()
{
  if (isRequesting) {
    return;
  }

  m2n->prepareEstablishment(localName, remoteName);
}

void BoundM2N::connectPrimaryRanks(std::string_view configHash)
{
  std::string fullLocalName = localName;

  if (isRequesting) {
    m2n->requestPrimaryRankConnection(remoteName, fullLocalName, configHash);
  } else {
    m2n->acceptPrimaryRankConnection(fullLocalName, remoteName, configHash);
  }
}

void BoundM2N::connectSecondaryRanks()
{
  if (m2n->usesTwoLevelInitialization()) {
    PRECICE_DEBUG("Update secondary connections");
    m2n->completeSecondaryRanksConnection();
  } else {
    if (isRequesting) {
      PRECICE_DEBUG("Awaiting secondary connections from {}", remoteName);
      m2n->requestSecondaryRanksConnection(remoteName, localName);
      PRECICE_DEBUG("Established secondary connections from {}", remoteName);
    } else {
      PRECICE_DEBUG("Establishing secondary connections to {}", remoteName);
      m2n->acceptSecondaryRanksConnection(localName, remoteName);
      PRECICE_DEBUG("Established  secondary connections to {}", remoteName);
    }
  }
}

void BoundM2N::preConnectSecondaryRanks()
{
  if (not m2n->usesTwoLevelInitialization())
    return;

  PRECICE_WARN("Two-level initialization is still in beta testing. Several edge cases are known to fail. Please report problems nevertheless.");
  Event e("bound-m2n.preConnectSecondaryRanks");

  // Accepting side (set up, gather connection info, communicate to requesting side)
  if (!isRequesting) {
    // Set up accepting side
    PRECICE_DEBUG("Setting up preliminary secondary connections from {}", localName);
    std::string connectionInfo;
    connectionInfo = m2n->prepareAcceptSecondaryRanksPreConnection(localName, remoteName);
    PRECICE_DEBUG("Set up preliminary secondary connections from {}. Ready for establishing connections.", localName);

    // Gather connection info and communicate it
    if (utils::IntraComm::isSecondary()) {
      Event e1("bound-m2n.gatherSendConnectionInfo");
      utils::IntraComm::getCommunication()->gatherRanges(connectionInfo);
      e1.stop();
    } else { // Primary
      // Gather connection info
      Event e1("bound-m2n.gatherConnectionInfoMap");

      std::map<Rank, std::string> connectionInfoMap;

      std::vector<std::string> allConnectionInfos = utils::IntraComm::getCommunication()->gatherRanges(connectionInfo);

      // Store the primary rank's connection info as well
      connectionInfoMap.emplace(0, allConnectionInfos[0]);

      for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
        PRECICE_ASSERT(secondaryRank < allConnectionInfos.size(), "Connection information from secondary rank not received");
        connectionInfoMap.emplace(secondaryRank, allConnectionInfos[secondaryRank]);
      }

      e1.stop();

      // Communicate connection info via primaries
      Event e2("bound-m2n.sendConnectionInfoMap");
      com::sendConnectionInfoMap(*m2n->getPrimaryRankCommunication(), 0, connectionInfoMap);
      e2.stop();
    }

    PRECICE_DEBUG("Establishing preliminary secondary connections to {}", remoteName);
    m2n->finishAcceptSecondaryRanksPreConnection(localName, remoteName);
    PRECICE_DEBUG("Established preliminary secondary connections to {}", remoteName);
  } else { // isRequesting
    com::serialize::SerializedConnectionInfoMap::ConnectionInfoMap connectionInfoMap;

    if (utils::IntraComm::isPrimary()) {
      // Communicate connection info via primaries
      Event e1("bound-m2n.receiveConnectionInfoMap");
      com::receiveConnectionInfoMap(*m2n->getPrimaryRankCommunication(), 0, connectionInfoMap);
      e1.stop();

      // Scatter connection info
      Event e2("bound-m2n.scatterSendConnectionInfoMap");
      com::broadcastSendConnectionInfoMap(*utils::IntraComm::getCommunication(), connectionInfoMap);
      e2.stop();
    } else {
      // Receive connection info from scatter
      Event e1("bound-m2n.scatterReceiveConnectionInfoMap");
      com::broadcastReceiveConnectionInfoMap(*utils::IntraComm::getCommunication(), connectionInfoMap);
      e1.stop();
    }

    // Connect to the accepting side
    PRECICE_DEBUG("Awaiting preliminary secondary connections from {}", remoteName);
    m2n->requestSecondaryRanksPreConnection(remoteName, localName, connectionInfoMap);
    PRECICE_DEBUG("Established preliminary secondary connections from {}", remoteName);
  }
}

void BoundM2N::cleanupEstablishment()
{
  if (isRequesting) {
    return;
  }
  waitForSecondaryRanks();
  if (!utils::IntraComm::isSecondary()) {
    m2n->cleanupEstablishment(localName, remoteName);
  }
}

void BoundM2N::waitForSecondaryRanks()
{
  if (utils::IntraComm::isPrimary()) {
    for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
      int item = 0;
      utils::IntraComm::getCommunication()->receive(item, rank);
      PRECICE_ASSERT(item > 0);
    }
  }
  if (utils::IntraComm::isSecondary()) {
    int item = utils::IntraComm::getRank();
    utils::IntraComm::getCommunication()->send(item, 0);
  }
}

} // namespace precice::m2n
