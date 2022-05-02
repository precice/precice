#include "GatherScatterCommunication.hpp"

#include "GatherScatterComFactory.hpp"

#include <utility>

namespace precice {
namespace m2n {
GatherScatterComFactory::GatherScatterComFactory(com::PtrCommunication intraComm) : _intraComm(std::move(intraComm)) {}

DistributedCommunication::SharedPointer GatherScatterComFactory::newDistributedCommunication(mesh::PtrMesh mesh)
{
  return DistributedCommunication::SharedPointer(new GatherScatterCommunication(_intraComm, mesh));
}
} // namespace m2n
} // namespace precice
