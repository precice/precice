#ifndef PRECICE_NO_MPI

#include "MPISinglePortsCommunicationFactory.hpp"
#include <memory>
#include <utility>

#include "MPISinglePortsCommunication.hpp"
#include "com/SharedPointer.hpp"

namespace precice::com {
MPISinglePortsCommunicationFactory::MPISinglePortsCommunicationFactory(std::string addressDirectory)
    : _addressDirectory(std::move(addressDirectory))
{
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

PtrCommunication MPISinglePortsCommunicationFactory::newCommunication()
{
  return std::make_shared<MPISinglePortsCommunication>(_addressDirectory);
}

std::string MPISinglePortsCommunicationFactory::addressDirectory()
{
  return _addressDirectory;
}
} // namespace precice::com

#endif // not PRECICE_NO_MPI
