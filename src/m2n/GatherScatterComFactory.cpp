#include "GatherScatterCommunication.hpp"

#include "GatherScatterComFactory.hpp"

#include <utility>

namespace precice {
namespace m2n {
GatherScatterComFactory::GatherScatterComFactory(
    com::PtrCommunication masterCom)
    : _masterCom(std::move(masterCom))
{
}

DistributedCommunication::SharedPointer
GatherScatterComFactory::newDistributedCommunication(mesh::PtrMesh mesh)
{
  return DistributedCommunication::SharedPointer(
      new GatherScatterCommunication(_masterCom, mesh));
}
} // namespace m2n
} // namespace precice
