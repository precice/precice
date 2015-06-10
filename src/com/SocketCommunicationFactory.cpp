// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_SOCKETS

#include "SocketCommunicationFactory.hpp"

#include "SocketCommunication.hpp"

namespace precice {
namespace com {
SocketCommunicationFactory::SocketCommunicationFactory(
    unsigned short portNumber,
    bool reuseAddress,
    std::string const& networkName,
    std::string const& addressDirectory)
    : _portNumber(portNumber)
    , _reuseAddress(reuseAddress)
    , _networkName(networkName)
    , _addressDirectory(addressDirectory) {
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

SocketCommunicationFactory::SocketCommunicationFactory(
    std::string const& addressDirectory)
    : SocketCommunicationFactory(0, false, "lo", addressDirectory) {
}

Communication::SharedPointer
SocketCommunicationFactory::newCommunication() {
  return Communication::SharedPointer(new SocketCommunication(
      _portNumber, _reuseAddress, _networkName, _addressDirectory));
}

std::string
SocketCommunicationFactory::addressDirectory() {
  return _addressDirectory;
}
}
} // namespace precice, com

#endif // not PRECICE_NO_SOCKETS
