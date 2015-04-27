// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#pragma once

#include "DistributedComFactory.hpp"

#include "com/CommunicationFactory.hpp"

namespace precice {
namespace m2n {
class PointToPointComFactory : public DistributedComFactory {
public:
  /**
   * @brief Constructor.
   */
  PointToPointComFactory(com::CommunicationFactory::SharedPointer comFactory);

  DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh);

private:
  // @brief communication factory for 1:M communications
  com::CommunicationFactory::SharedPointer _comFactory;
};
}
} // namespace precice, m2n
