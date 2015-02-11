// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#pragma once

#include "SharedPointer.hpp"

namespace precice {
namespace m2n {
class DistributedCommunicationFactory {
public:
  /**
   * @brief Destructor.
   */
  virtual ~DistributedCommunicationFactory() {};

  virtual PtrDistributedCommunication newDistributedCommunication(mesh::PtrMesh mesh) = 0;
};

}} // namespace precice, m2n

