#include "partition/ProvidedPartition.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/M2N.hpp"

namespace precice {
namespace partition {

logging::Logger ProvidedPartition:: _log ( "precice::partition::ProvidedPartition" );

ProvidedPartition::ProvidedPartition
(
    bool hasToSend)
:
    Partition (),
    _hasToSend(hasToSend)
{}

void ProvidedPartition::communicate()
{
  TRACE();

  if(_hasToSend){

    // Temporary globalMesh such that the master also keeps his local mesh
    mesh::Mesh globalMesh(_mesh->getName(), _mesh->getDimensions(), _mesh->isFlipNormals());

    if( not utils::MasterSlave::_slaveMode ){
      globalMesh.addMesh(*_mesh); //add local master mesh to global mesh
    }

    // Gather Mesh
    INFO("Gather mesh " << _mesh->getName() );
    if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode ) {
      if (utils::MasterSlave::_slaveMode) {
        com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh( *_mesh, 0 );
      }
      else{ // Master
        assertion(utils::MasterSlave::_rank==0);
        assertion(utils::MasterSlave::_size>1);

        for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
          com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh ( globalMesh, rankSlave);
          DEBUG("Received sub-mesh, from slave: " << rankSlave <<", global vertexCount: " << globalMesh.vertices().size());
        }
      }
    }

    //set global index
    if (not utils::MasterSlave::_slaveMode) {
       for(mesh::Vertex& v : globalMesh.vertices()){
         v.setGlobalIndex(v.getID());
       }
    }

    // Send (global) Mesh
    INFO("Send global mesh " << _mesh->getName());
    if (not utils::MasterSlave::_slaveMode) {
      CHECK ( globalMesh.vertices().size() > 0, "The provided mesh " << globalMesh.getName() << " is invalid (possibly empty).");
      com::CommunicateMesh(_m2n->getMasterCommunication()).sendMesh ( globalMesh, 0 );
    }

  } //_hasToSend
}

void ProvidedPartition::compute()
{
  TRACE();

  int numberOfVertices = _mesh->vertices().size();

  if (utils::MasterSlave::_slaveMode) {
    utils::MasterSlave::_communication->send(numberOfVertices,0);
  }
  else if(utils::MasterSlave::_masterMode) {
    assertion(utils::MasterSlave::_size>1);
    int vertexCounter = 0;

    //add master vertices
    for(int i=0; i<numberOfVertices; i++){
      _mesh->getVertexDistribution()[0].push_back(vertexCounter);
      vertexCounter++;
    }

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      utils::MasterSlave::_communication->receive(numberOfVertices,rankSlave);
      for(int i=0; i<numberOfVertices; i++){
        _mesh->getVertexDistribution()[rankSlave].push_back(vertexCounter);
        vertexCounter++;
      }
    }
    _mesh->setGlobalNumberOfVertices(vertexCounter);
  }
  else{ //coupling mode
    _mesh->setGlobalNumberOfVertices(_mesh->vertices().size());
  }


  // TODO kann man computeDistr nicht direkt hier machen?
  // TODO owner berechnen
  // TODO globalIndex an slaves schicken
  // TODO globalIndex setzen falls nicht gesendet wurde
  _mesh->computeDistribution();
}

}} // namespace precice, partition
