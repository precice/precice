#include "partition/Partition.hpp"
#include "com/Communication.hpp"
#include "mesh/Mesh.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace partition {

Partition::Partition(mesh::PtrMesh mesh)
    : _mesh(mesh)
{
}

void Partition::computeVertexOffsets()
{
  PRECICE_TRACE();
  PRECICE_DEBUG("Generate vertex offsets");
  if (utils::MasterSlave::isSlave()) {
    utils::MasterSlave::_communication->broadcast(_mesh->getVertexOffsets(), 0);
    PRECICE_DEBUG("My vertex offsets: " << _mesh->getVertexOffsets());

  } else if (utils::MasterSlave::isMaster()) {
    _mesh->getVertexOffsets().resize(utils::MasterSlave::getSize());
    _mesh->getVertexOffsets()[0] = _mesh->getVertexDistribution()[0].size();
    for (int rank = 1; rank < utils::MasterSlave::getSize(); rank++) {
      _mesh->getVertexOffsets()[rank] = _mesh->getVertexDistribution()[rank].size() + _mesh->getVertexOffsets()[rank - 1];
    }
    PRECICE_DEBUG("My vertex offsets: " << _mesh->getVertexOffsets());
    utils::MasterSlave::_communication->broadcast(_mesh->getVertexOffsets());

  } else { //coupling mode
    _mesh->getVertexOffsets().push_back(_mesh->vertices().size());
  }
}

} // namespace partition
} // namespace precice
