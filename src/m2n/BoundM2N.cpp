#include <memory>

#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/BoundM2N.hpp"
#include "m2n/M2N.hpp"
#include "precice/types.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace m2n {

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
    m2n->requestPrimaryConnection(remoteName, fullLocalName);
  } else {
    m2n->acceptPrimaryConnection(fullLocalName, remoteName);
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
      m2n->requestSecondaryConnections(remoteName, localName);
      PRECICE_DEBUG("Established secondary connections from {}", remoteName);
    } else {
      PRECICE_DEBUG("Establishing secondary connections to {}", remoteName);
      m2n->acceptSecondaryConnections(localName, remoteName);
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
    m2n->acceptSecondaryPreConnections(localName, remoteName);
    PRECICE_DEBUG("Established preliminary secondary connections to {}", remoteName);
  }
}

void BoundM2N::cleanupEstablishment()
{
  if (isRequesting) {
    return;
  }
  waitForSecondaryRanks();
  if (!utils::MasterSlave::isSecondary()) {
    m2n->cleanupEstablishment(localName, remoteName);
  }
}

void BoundM2N::waitForSecondaryRanks()
{
  if (utils::MasterSlave::isPrimary()) {
    for (Rank rank : utils::MasterSlave::allSecondaryRanks()) {
      int item = 0;
      utils::MasterSlave::getCommunication()->receive(item, rank);
      PRECICE_ASSERT(item > 0);
    }
  }
  if (utils::MasterSlave::isSecondary()) {
    int item = utils::MasterSlave::getRank();
    utils::MasterSlave::getCommunication()->send(item, 0);
  }
}

} // namespace m2n
} // namespace precice
