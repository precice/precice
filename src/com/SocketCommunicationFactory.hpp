// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_SOCKET_COMMUNICATION_FACTORY_HPP_
#define PRECICE_COM_SOCKET_COMMUNICATION_FACTORY_HPP_

#include "CommunicationFactory.hpp"

#include <string>

namespace precice {
namespace com {
class SocketCommunicationFactory : public CommunicationFactory {
public:
  /**
   * @brief Constructor.
   */
  SocketCommunicationFactory(
      std::string const& networkName,
      unsigned short portNumber,
      std::string const& ipAddressExchangeDirectory = "");

  /**
   * @brief Constructor.
   */
  SocketCommunicationFactory(
      std::string const& ipAddressExchangeDirectory = "");

  PtrCommunication newCommunication();

private:
  std::string _networkName;
  unsigned short _portNumber;
  std::string _ipAddressExchangeDirectory;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_SOCKET_COMMUNICATION_FACTORY_HPP_ */
