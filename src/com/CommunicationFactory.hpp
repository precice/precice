#pragma once

#include <stdexcept>
#include "com/SharedPointer.hpp"

namespace precice
{
namespace com
{
class CommunicationFactory
{

public:
  virtual ~CommunicationFactory(){};

  virtual PtrCommunication newCommunication() = 0;

  virtual std::string addressDirectory()
  {
    throw std::runtime_error("Not available!");
  }
};
}
} // namespace precice, com
