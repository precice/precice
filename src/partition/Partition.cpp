#include "partition/Partition.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Globals.hpp"
#include "com/Communication.hpp"
#include "utils/MasterSlave.hpp"


namespace precice {
namespace partition {

logging::Logger Partition::_log("partition::Partition" );

Partition:: Partition
(mesh::PtrMesh mesh)
:
  _mesh(mesh),
  _fromMapping(),
  _toMapping(),
  _m2n()
{}


void Partition:: computeVertexOffsets(){
  TRACE();
  DEBUG("Generate vertex offsets");
  if (utils::MasterSlave::_slaveMode) {
    _mesh->getVertexOffsets().resize(utils::MasterSlave::_size);
    utils::MasterSlave::_communication->broadcast(_mesh->getVertexOffsets().data(),_mesh->getVertexOffsets().size(),0);
    DEBUG("My vertex offsets: " << _mesh->getVertexOffsets());
  }
  else if (utils::MasterSlave::_masterMode) {
    _mesh->getVertexOffsets().resize(utils::MasterSlave::_size);
    _mesh->getVertexOffsets()[0] = _mesh->getVertexDistribution()[0].size();
    for (int rank = 1; rank < utils::MasterSlave::_size; rank++){
      _mesh->getVertexOffsets()[rank] = _mesh->getVertexDistribution()[rank].size() + _mesh->getVertexOffsets()[rank-1];
    }
    DEBUG("My vertex offsets: " << _mesh->getVertexOffsets());
    utils::MasterSlave::_communication->broadcast(_mesh->getVertexOffsets().data(),_mesh->getVertexOffsets().size());
  }
  else{ //coupling mode
    _mesh->getVertexOffsets().push_back(_mesh->vertices().size());
  }
}


}}
