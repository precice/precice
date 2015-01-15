// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_COMMUNICATION_FACTORY_HPP_
#define PRECICE_COM_COMMUNICATION_FACTORY_HPP_

#include "SharedPointer.hpp"

namespace precice {
namespace com {
class CommunicationFactory {
public:
  /**
   * @brief Destructor.
   */
  virtual ~CommunicationFactory() = 0;

  virtual PtrCommunication newCommunication() = 0;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_COMMUNICATION_FACTORY_HPP_ */
