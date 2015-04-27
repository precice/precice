// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "GatherScatterCommunication.hpp"

#include "GatherScatterComFactory.hpp"

namespace precice {
namespace m2n {
GatherScatterComFactory::GatherScatterComFactory(
    com::Communication::SharedPointer masterCom)
    : _masterCom(masterCom) {
}

DistributedCommunication::SharedPointer
GatherScatterComFactory::newDistributedCommunication(mesh::PtrMesh mesh) {
  return DistributedCommunication::SharedPointer(
      new GatherScatterCommunication(_masterCom, mesh));
}
}
} // namespace precice, m2n
