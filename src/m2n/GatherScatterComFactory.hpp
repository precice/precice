#pragma once

#include "DistributedComFactory.hpp"
#include "com/SharedPointer.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice::m2n {
class GatherScatterComFactory : public DistributedComFactory {
public:
  GatherScatterComFactory(com::PtrCommunication intraComm);

  DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh) override;

private:
  /// communication between the primary processes
  com::PtrCommunication _intraComm;
};
} // namespace precice::m2n
