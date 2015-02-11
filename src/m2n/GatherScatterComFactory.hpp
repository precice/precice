// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#pragma once

#include "DistributedCommunicationFactory.hpp"

namespace precice {
namespace m2n {
class GatherScatterComFactory : public DistributedCommunicationFactory {
public:
  /**
   * @brief Constructor.
   */
  GatherScatterComFactory(com::PtrCommunication masterCom);


  PtrDistributedCommunication newDistributedCommunication(mesh::PtrMesh mesh);

private:

  // @brief communication between the master processes
  com::PtrCommunication _masterCom;
};

}} // namespace precice, m2n

