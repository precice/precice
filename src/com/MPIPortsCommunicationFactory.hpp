#ifndef PRECICE_NO_MPI

#pragma once

#include "CommunicationFactory.hpp"
#include "com/SharedPointer.hpp"

#include <string>

namespace precice
{
namespace com
{
class MPIPortsCommunicationFactory : public CommunicationFactory
{
public:
  explicit MPIPortsCommunicationFactory(std::string const &addressDirectory = ".");

  PtrCommunication newCommunication();

  std::string addressDirectory();

private:
  std::string _addressDirectory;
};
}
} // namespace precice, com

#endif // not PRECICE_NO_MPI
