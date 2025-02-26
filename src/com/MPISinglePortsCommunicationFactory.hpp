#ifndef PRECICE_NO_MPI

#pragma once

#include "CommunicationFactory.hpp"
#include "com/SharedPointer.hpp"

#include <string>

namespace precice::com {
class MPISinglePortsCommunicationFactory : public CommunicationFactory {
public:
  explicit MPISinglePortsCommunicationFactory(std::string addressDirectory = ".");

  PtrCommunication newCommunication() override;

  std::string addressDirectory() override;

private:
  std::string _addressDirectory;
};
} // namespace precice::com

#endif // not PRECICE_NO_MPI
