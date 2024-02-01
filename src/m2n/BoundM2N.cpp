#include <memory>

#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/BoundM2N.hpp"
#include "m2n/M2N.hpp"
#include "precice/impl/Types.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::m2n {

void BoundM2N::prepareEstablishment()
{
  if (isRequesting) {
    return;
  }

  m2n->prepareEstablishment(localName, remoteName);
}

void BoundM2N::connectPrimaryRanks()
{
  std::string fullLocalName = localName;

  if (isRequesting) {
    m2n->requestPrimaryRankConnection(remoteName, fullLocalName);
  } else {
    m2n->acceptPrimaryRankConnection(fullLocalName, remoteName);
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

  if (isRequesting) {
    PRECICE_DEBUG("Awaiting preliminary secondary connections from {}", remoteName);
    m2n->requestSecondaryRanksPreConnection(remoteName, localName);
    PRECICE_DEBUG("Established preliminary secondary connections from {}", remoteName);
  } else {
    PRECICE_DEBUG("Establishing preliminary secondary connections to {}", remoteName);
    m2n->acceptSecondaryRanksPreConnection(localName, remoteName);
    PRECICE_DEBUG("Established preliminary secondary connections to {}", remoteName);
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
