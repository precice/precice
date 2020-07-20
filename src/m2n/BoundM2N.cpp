#include "m2n/BoundM2N.hpp"
#include <memory>
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/M2N.hpp"
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

void BoundM2N::connectMasters()
{
  std::string fullLocalName = localName;

  if (isRequesting) {
    m2n->requestMasterConnection(remoteName, fullLocalName);
  } else {
    m2n->acceptMasterConnection(fullLocalName, remoteName);
  }
}

void BoundM2N::connectSlaves()
{
  if (m2n->usesTwoLevelInitialization()) {
    PRECICE_DEBUG("Update slaves connections");
    m2n->completeSlavesConnection();
  } else {
    if (isRequesting) {
      PRECICE_DEBUG("Awaiting slaves connection from " << remoteName);
      m2n->requestSlavesConnection(remoteName, localName);
      PRECICE_DEBUG("Established slaves connection from " << remoteName);
    } else {
      PRECICE_DEBUG("Establishing slaves connection to " << remoteName);
      m2n->acceptSlavesConnection(localName, remoteName);
      PRECICE_DEBUG("Established  slaves connection to " << remoteName);
    }
  }
}

void BoundM2N::preConnectSlaves()
{
  if (not m2n->usesTwoLevelInitialization())
    return;

  PRECICE_WARN("Two-level initialization is still in beta testing. Several edge cases are known to fail. Please report problems nevertheless.");

  if (isRequesting) {
    PRECICE_DEBUG("Awaiting preliminary slaves connection from " << remoteName);
    m2n->requestSlavesPreConnection(remoteName, localName);
    PRECICE_DEBUG("Established preliminary slaves connection from " << remoteName);
  } else {
    PRECICE_DEBUG("Establishing preliminary slaves connection to " << remoteName);
    m2n->acceptSlavesPreConnection(localName, remoteName);
    PRECICE_DEBUG("Established preliminary slaves connection to " << remoteName);
  }
}

void BoundM2N::cleanupEstablishment()
{
  if (isRequesting) {
    return;
  }
  waitForSlaves();
  if (!utils::MasterSlave::isSlave()) {
    m2n->cleanupEstablishment(localName, remoteName);
  }
}

void BoundM2N::waitForSlaves()
{
  if (utils::MasterSlave::isMaster()) {
    for (int rank = 1; rank < utils::MasterSlave::getSize(); ++rank) {
      int item = 0;
      utils::MasterSlave::_communication->receive(item, rank);
      PRECICE_ASSERT(item > 0);
    }
  }
  if (utils::MasterSlave::isSlave()) {
    int item = utils::MasterSlave::getRank();
    utils::MasterSlave::_communication->send(item, 0);
  }
}

} // namespace m2n
} // namespace precice
