#include "partition/ProvidedBoundingBox.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EventUtils.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

using precice::utils::Event;

namespace precice
{
namespace partition
{

ProvidedBoundingBox::ProvidedBoundingBox(mesh::PtrMesh mesh,
                                         bool          hasToSend,
                                         double        safetyFactor)
    : Partition(mesh),
      _hasToSend(hasToSend),
      _dimensions(mesh->getDimensions()),
      _safetyFactor(safetyFactor)
{
}

void ProvidedBoundingBox::communicateBoundingBox()
{
  PRECICE_TRACE();

  if (!_hasToSend)
    return;

  // each rank sends its bb to master
  if (utils::MasterSlave::isSlave()) { //slave
    com::CommunicateBoundingBox(utils::MasterSlave::_communication).sendBoundingBox(_mesh->getBoundingBox(), 0);
  } else { // Master

    PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
    PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

    // to store the collection of bounding boxes
    mesh::Mesh::BoundingBoxMap bbm;

    // master stores its bb into bbm
    bbm[0] = _mesh->getBoundingBox();

    // master receives bbs from slaves and stores them in bbm
    for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
      // initialize bbm
      bbm[rankSlave] = mesh::Mesh::BoundingBox(_dimensions);
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).receiveBoundingBox(bbm[rankSlave], rankSlave);
    }

    // master sends number of ranks and bbm to the other master
    _m2ns[0]->getMasterCommunication()->send(utils::MasterSlave::getSize(), 0);
    com::CommunicateBoundingBox(_m2ns[0]->getMasterCommunication()).sendBoundingBoxMap(bbm, 0);
  }
}

void ProvidedBoundingBox::computeBoundingBox()
{
  if (!_hasToSend)
    return;
  
  PRECICE_TRACE();

  // size of the feedbackmap
  int remoteConnectionMapSize = 0;
  std::vector<int> connectedRanksList;

  std::map<int, std::vector<int>> remoteConnectionMap;

  if (utils::MasterSlave::isMaster()) { 
    PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

    // master receives feedback map (map of other participant ranks -> connected ranks at this participant)
    // from other participants master
    _m2ns[0]->getMasterCommunication()->receive(connectedRanksList, 0);
    remoteConnectionMapSize = connectedRanksList.size();
    
    for (auto &rank : connectedRanksList) {
      remoteConnectionMap[rank] = {-1};
    }
    if (remoteConnectionMapSize != 0)
      com::CommunicateBoundingBox(_m2ns[0]->getMasterCommunication()).receiveConnectionMap(remoteConnectionMap, 0);

    // broadcast the received feedbackMap
    utils::MasterSlave::_communication->broadcast(connectedRanksList);
    if (remoteConnectionMapSize != 0) {
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendConnectionMap(remoteConnectionMap);
    }

    // master checks which ranks are connected to it
    for (auto &remoteRank : remoteConnectionMap) {
      for (auto &includedRank : remoteRank.second) {
        if (utils::MasterSlave::getRank() == includedRank) {
          _mesh->getConnectedRanks().push_back(remoteRank.first);
        }
      }
    }

  } else { // Slave

    utils::MasterSlave::_communication->broadcast(connectedRanksList, 0);

    if (!connectedRanksList.empty())
    {
      for (auto &rank : connectedRanksList) {
        remoteConnectionMap[rank] = {-1};
      }
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveConnectionMap(remoteConnectionMap);
    }

    for (auto &remoteRank : remoteConnectionMap) {
      for (auto &includedRanks : remoteRank.second) {
        if (utils::MasterSlave::getRank() == includedRanks) {
          _mesh->getConnectedRanks().push_back(remoteRank.first);
        }
      }
    }
  }
}

void ProvidedBoundingBox::communicate()
{

  if (!_hasToSend)
    return;  
  /*
   * First we should set global index for each vertex
   * This global vertex id is needed for final filtering in 
   * the received partition
   */

  // the maximum vertex global index for each rank
  int vertexMaxGlobalID = 0;
  // the minimum vertex global index for each rank
  int vertexMinGlobalID = 0;
  
  // offset to set global vertex id for each rank vertices
  int offset = 0;

  if (utils::MasterSlave::isMaster()) {//Master
    int vertexCounter = 0;

    // set global indexes for master rank mesh partition 
    for(int i=0; i<(int)_mesh->vertices().size(); i++){
      _mesh->getVertexDistribution()[0].push_back(vertexCounter);
      _mesh->vertices()[i].setGlobalIndex(vertexCounter);
      vertexCounter++;
    }

    // in master rank max global id is equal to number of vertices-1 
    vertexMaxGlobalID = _mesh->vertices().size()-1;
    
    // receive number of vertices for each rank at master to produce offset for the same rank
    // @todo: This is a gather operation which should be replaced with a better solution 
    int numberOfVertices = -1;
    for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++)
    {
      utils::MasterSlave::_communication->receive(numberOfVertices,rankSlave);
      utils::MasterSlave::_communication->send(vertexCounter,rankSlave);
      for (int i = 0; i < numberOfVertices; i++) {
        _mesh->getVertexDistribution()[rankSlave].push_back(vertexCounter);
        vertexCounter++;
      }
    }

    // set the total number of vertices for this participant
    _mesh->setGlobalNumberOfVertices(vertexCounter);
    // broadcast the total number of vertices to slaves
    utils::MasterSlave::_communication->broadcast(vertexCounter);
    // communicate the total number of vertices to the other participants master 
    _m2ns[0]->getMasterCommunication()->send(vertexCounter, 0);
  }
  else
  {
    // send number of vertices to master
    utils::MasterSlave::_communication->send((int)_mesh->vertices().size() ,0);
    
    // receive the offset from master
    utils::MasterSlave::_communication->receive(offset ,0);    
    vertexMinGlobalID = offset;

    // set global indexes for each vertex
    for(int i=0; i<(int)_mesh->vertices().size(); i++)
    {
      _mesh->vertices()[i].setGlobalIndex(offset+i);
    }           
     
    // set the max global index
    vertexMaxGlobalID = offset + _mesh->vertices().size() -1;

    // receive global number of vertices for this participant from master
    int globalNumberOfVertices = -1;
    utils::MasterSlave::_communication->broadcast(globalNumberOfVertices, 0);
    PRECICE_ASSERT(globalNumberOfVertices != -1);
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);   
  } 

  // each rank sends its min/max global vertex index to connected remote ranks
  _m2ns[0]->broadcastSend(vertexMinGlobalID, *_mesh);
  _m2ns[0]->broadcastSend(vertexMaxGlobalID, *_mesh);

  // each rank sends its mesh partition to connected remote ranks
  _m2ns[0]->broadcastSendLocalMesh(*_mesh);
  
  createOwnerInformation();  
  computeVertexOffsets();
}

void ProvidedBoundingBox::compute()
{
  if (!_hasToSend)
    return;

  // receive communication map from all remote connected ranks
   _m2ns[0]->broadcastReceiveLCM(_mesh->getCommunicationMap(), *_mesh);    
}

void ProvidedBoundingBox::createOwnerInformation()
{
  PRECICE_TRACE();
  for (mesh::Vertex &v : _mesh->vertices()) {
    v.setOwner(true);
  }
}

} // namespace partition
} // namespace precice
