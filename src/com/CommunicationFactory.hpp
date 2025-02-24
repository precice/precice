#pragma once

#include <stdexcept>
#include <string>
#include "com/SharedPointer.hpp"

namespace precice::com {
class CommunicationFactory {

public:
  virtual ~CommunicationFactory() = default;

  virtual PtrCommunication newCommunication() = 0;

  virtual std::string addressDirectory()
  {
    throw std::runtime_error("Not available!");
  }
};
} // namespace precice::com
