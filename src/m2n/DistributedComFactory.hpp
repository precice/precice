// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#pragma once

#include "DistributedCommunication.hpp"

#include <memory>

namespace precice {
namespace m2n {
class DistributedComFactory {
public:
  using SharedPointer = std::shared_ptr<DistributedComFactory>;

public:
  /**
   * @brief Destructor.
   */
  virtual ~DistributedComFactory(){};

  virtual DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh) = 0;
};
}
} // namespace precice, m2n
