#ifndef PRECICE_NO_MPI

#include "MPISinglePortsCommunicationFactory.hpp"
#include <memory>
#include "MPISinglePortsCommunication.hpp"
#include "com/SharedPointer.hpp"

namespace precice {
namespace com {
MPISinglePortsCommunicationFactory::MPISinglePortsCommunicationFactory(std::string const &addressDirectory)
    : _addressDirectory(addressDirectory)
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
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
