#ifndef PRECICE_NO_SOCKETS

#pragma once

#include "CommunicationFactory.hpp"
#include "com/SharedPointer.hpp"

#include <string>

namespace precice
{
namespace com
{
class SocketCommunicationFactory : public CommunicationFactory
{
public:
  SocketCommunicationFactory(unsigned short     portNumber       = 0,
                             bool               reuseAddress     = false,
                             std::string const &networkName      = "lo",
                             std::string const &addressDirectory = ".");

  explicit SocketCommunicationFactory(std::string const &addressDirectory);

  PtrCommunication newCommunication();

  std::string addressDirectory();

private:
  unsigned short _portNumber;
  bool           _reuseAddress;
  std::string    _networkName;
  std::string    _addressDirectory;
};
}
} // namespace precice, com

#endif // not PRECICE_NO_SOCKETS
