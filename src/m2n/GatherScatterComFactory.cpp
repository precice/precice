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
