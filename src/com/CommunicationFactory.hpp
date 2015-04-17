// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_COMMUNICATION_FACTORY_HPP_
#define PRECICE_COM_COMMUNICATION_FACTORY_HPP_

#include "Communication.hpp"

#include <memory>
#include <stdexcept>

namespace precice {
namespace com {
class CommunicationFactory {
public:
  using SharedPointer = std::shared_ptr<CommunicationFactory>;

public:
  /**
   * @brief Destructor.
   */
  virtual ~CommunicationFactory(){};

  virtual Communication::SharedPointer newCommunication() = 0;

  virtual std::string
  addressDirectory() {
    throw std::runtime_error("Not available!");
  }
};
}
} // namespace precice, com

#endif /* PRECICE_COM_COMMUNICATION_FACTORY_HPP_ */
