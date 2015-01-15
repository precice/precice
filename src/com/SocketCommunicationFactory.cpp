// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "SocketCommunication.hpp"

#include "SocketCommunicationFactory.hpp"

namespace precice {
namespace com {
SocketCommunicationFactory::SocketCommunicationFactory(
    std::string const& networkName,
    unsigned short portNumber,
    std::string const& ipAddressExchangeDirectory)
    : _networkName(networkName)
    , _portNumber(portNumber)
    , _ipAddressExchangeDirectory(ipAddressExchangeDirectory) {
}

SocketCommunicationFactory::SocketCommunicationFactory(
    std::string const& ipAddressExchangeDirectory)
    : _ipAddressExchangeDirectory(ipAddressExchangeDirectory) {
}

PtrCommunication
SocketCommunicationFactory::newCommunication() {
  return PtrCommunication(new SocketCommunication(
      _networkName, _portNumber, _ipAddressExchangeDirectory));
}
}
} // namespace precice, com
