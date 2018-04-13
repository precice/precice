#pragma once

#include "DistributedComFactory.hpp"

namespace precice
{
namespace m2n
{
class GatherScatterComFactory : public DistributedComFactory
{
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
