#pragma once

#include "CommunicationFactory.hpp"
#include "com/SharedPointer.hpp"
#include "utils/networking.hpp"

#include <string>

namespace precice::com {
class SocketCommunicationFactory : public CommunicationFactory {
public:
  SocketCommunicationFactory(unsigned short portNumber       = 0,
                             bool           reuseAddress     = false,
                             std::string    networkName      = utils::networking::loopbackInterfaceName(),
                             std::string    addressDirectory = ".");

  explicit SocketCommunicationFactory(std::string const &addressDirectory);

  PtrCommunication newCommunication() override;

  std::string addressDirectory() override;

private:
  unsigned short _portNumber;
  bool           _reuseAddress;
  std::string    _networkName;
  std::string    _addressDirectory;
};
} // namespace precice::com
