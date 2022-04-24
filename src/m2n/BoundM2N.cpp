#include <memory>

#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/BoundM2N.hpp"
#include "m2n/M2N.hpp"
#include "precice/types.hpp"
#include "utils/IntraComm.hpp"
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

void BoundM2N::connectPrimaries()
{
  std::string fullLocalName = localName;

  if (isRequesting) {
    m2n->requestPrimaryConnection(remoteName, fullLocalName);
  } else {
    m2n->acceptPrimaryConnection(fullLocalName, remoteName);
  }
}

void BoundM2N::connectSecondaries()
{
  if (m2n->usesTwoLevelInitialization()) {
    PRECICE_DEBUG("Update slaves connections");
    m2n->completeSecondariesConnection();
  } else {
    if (isRequesting) {
      PRECICE_DEBUG("Awaiting slaves connection from {}", remoteName);
      m2n->requestSecondariesConnection(remoteName, localName);
      PRECICE_DEBUG("Established slaves connection from {}", remoteName);
    } else {
      PRECICE_DEBUG("Establishing slaves connection to {}", remoteName);
      m2n->acceptSecondariesConnection(localName, remoteName);
      PRECICE_DEBUG("Established  slaves connection to {}", remoteName);
    }
  }
}

void BoundM2N::preConnectSecondaries()
{
  if (not m2n->usesTwoLevelInitialization())
    return;

  PRECICE_WARN("Two-level initialization is still in beta testing. Several edge cases are known to fail. Please report problems nevertheless.");

  if (isRequesting) {
    PRECICE_DEBUG("Awaiting preliminary slaves connection from {}", remoteName);
    m2n->requestSecondariesPreConnection(remoteName, localName);
    PRECICE_DEBUG("Established preliminary slaves connection from {}", remoteName);
  } else {
    PRECICE_DEBUG("Establishing preliminary slaves connection to {}", remoteName);
    m2n->acceptSecondariesPreConnection(localName, remoteName);
    PRECICE_DEBUG("Established preliminary slaves connection to {}", remoteName);
  }
}

void BoundM2N::cleanupEstablishment()
{
  if (isRequesting) {
    return;
  }
  waitForSecondaries();
  if (!utils::IntraComm::isSecondary()) {
    m2n->cleanupEstablishment(localName, remoteName);
  }
}

void BoundM2N::waitForSecondaries()
{
  if (utils::IntraComm::isPrimary()) {
    for (Rank rank : utils::IntraComm::allSecondaries()) {
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

} // namespace m2n
} // namespace precice
