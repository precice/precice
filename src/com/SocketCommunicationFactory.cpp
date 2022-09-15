#include "SocketCommunicationFactory.hpp"
#include <memory>
#include <utility>

#include "SocketCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "utils/networking.hpp"

namespace precice::com {
SocketCommunicationFactory::SocketCommunicationFactory(
    unsigned short portNumber,
    bool           reuseAddress,
    std::string    networkName,
    std::string    addressDirectory)
    : _portNumber(portNumber),
      _reuseAddress(reuseAddress),
      _networkName(std::move(networkName)),
      _addressDirectory(std::move(addressDirectory))
{
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

SocketCommunicationFactory::SocketCommunicationFactory(
    std::string const &addressDirectory)
    : SocketCommunicationFactory(0, false, utils::networking::loopbackInterfaceName(), addressDirectory)
{
}

PtrCommunication SocketCommunicationFactory::newCommunication()
{
  return std::make_shared<SocketCommunication>(
      _portNumber, _reuseAddress, _networkName, _addressDirectory);
}

std::string SocketCommunicationFactory::addressDirectory()
{
  return _addressDirectory;
}
} // namespace precice::com
