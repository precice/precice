#include "GatherScatterCommunication.hpp"

#include "GatherScatterComFactory.hpp"

#include <utility>

namespace precice {
namespace m2n {
GatherScatterComFactory::GatherScatterComFactory(
    com::PtrCommunication primaryCom)
    : _primaryCom(std::move(primaryCom))
{
}

DistributedCommunication::SharedPointer
GatherScatterComFactory::newDistributedCommunication(mesh::PtrMesh mesh)
{
  return DistributedCommunication::SharedPointer(
      new GatherScatterCommunication(_primaryCom, mesh));
}
} // namespace m2n
} // namespace precice
