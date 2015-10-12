// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_SOCKETS

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
  SocketCommunicationFactory(unsigned short portNumber = 0,
                             bool reuseAddress = false,
                             std::string const& networkName = "lo",
                             std::string const& addressDirectory = ".");

  /**
   * @brief Constructor.
   */
  SocketCommunicationFactory(std::string const& addressDirectory);

  Communication::SharedPointer newCommunication();

  std::string addressDirectory();

private:
  unsigned short _portNumber;
  bool _reuseAddress;
  std::string _networkName;
  std::string _addressDirectory;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_SOCKET_COMMUNICATION_FACTORY_HPP_ */

#endif // not PRECICE_NO_SOCKETS
