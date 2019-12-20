#ifndef PRECICE_NO_MPI

#pragma once

#include "CommunicationFactory.hpp"
#include "com/SharedPointer.hpp"

#include <string>

namespace precice {
namespace com {
class MPIPortsCommunicationFactory : public CommunicationFactory {
public:
  explicit MPIPortsCommunicationFactory(std::string const &addressDirectory = ".");

  PtrCommunication newCommunication() override;

  std::string addressDirectory() override;

private:
  std::string _addressDirectory;
};
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
