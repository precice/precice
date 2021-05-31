#pragma once

#include "DistributedComFactory.hpp"
#include "com/SharedPointer.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace m2n {
class GatherScatterComFactory : public DistributedComFactory {
public:
  GatherScatterComFactory(com::PtrCommunication masterCom);

  DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh);

private:
  /// communication between the master processes
  com::PtrCommunication _masterCom;
};
} // namespace m2n
} // namespace precice
