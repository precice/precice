#include "partition/ProvidedBoundingBox.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/M2N.hpp"
#include "utils/EventTimings.hpp"
#include "utils/Parallel.hpp"
#include "mesh/Mesh.hpp"
#include "mapping/Mapping.hpp"

using precice::utils::Event;

namespace precice {
namespace partition {

ProvidedBoundingBox::ProvidedBoundingBox
(mesh::PtrMesh mesh,
 bool hasToSend,
 double safetyFactor)
:
    Partition (mesh),
    _hasToSend(hasToSend),
    _dimensions(mesh->getDimensions()),
    _safetyFactor(safetyFactor)
{}

void ProvidedBoundingBox::communicateBoundingBox()
{

  if (_hasToSend) {

    mesh::Mesh::BoundingBox boundingBox; // use this to gather the bbs in the master rank
    // each rank sends its bb to master
    if (utils::MasterSlave::_slaveMode) {//slave
      boundingBox = _mesh->getBoundingBox();
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).sendBoundingBox(boundingBox, 0); 
    }
    else
    { // Master
      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);

      // master stores its bb into gloabalBB
      mesh::Mesh::BoundingBoxMap globalBB;
      globalBB[0] = _mesh->getBoundingBox();

      // master receives bbs from slaves and store them in global bb
      if (utils::MasterSlave::_size>1) {  
        for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
          com::CommunicateBoundingBox(utils::MasterSlave::_communication).receiveBoundingBox(boundingBox, rankSlave);       
          DEBUG("From slave " << rankSlave << ", bounding mesh: " << boundingBox[0].first
                << ", " << boundingBox[0].second << " and " << boundingBox[1].first << ", " << boundingBox[1].second);        
          globalBB[rankSlave] = boundingBox;
        }
      }
      
      //master sends set of boundingboxes (globalBB) to the other master
      _m2n->getMasterCommunication()->send(utils::MasterSlave::_size , 0);
      com::CommunicateBoundingBox(_m2n->getMasterCommunication()).sendBoundingBoxMap(globalBB,0);
    }                  
  }   
}

void ProvidedBoundingBox::computeBoundingBox()
{
  int numberOfVertices = _mesh->vertices().size();
  int remoteParComSize = -1;
  mesh::Mesh::FeedbackMap receivedFeedbackMap;
  if (not utils::MasterSlave::_slaveMode) {//Master
    assertion(utils::MasterSlave::_size>1);
    //master receives other partition communicator size and also a feedback which is a map:  list of other participant ranks -> connected ranks at this participant  
    _m2n->getMasterCommunication()->receive(remoteParComSize, 0);
    utils::MasterSlave::_communication->broadcast(remoteParComSize);
    for (int i=0; i < remoteParComSize; i++)
    {
      std::vector<int> initialFeedback;
      initialFeedback.push_back(-1);
      receivedFeedbackMap[i]=initialFeedback;
    }
    com::CommunicateBoundingBox(_m2n->getMasterCommunication()).receiveFeedbackMap(receivedFeedbackMap, 0 ); 
    com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendFeedbackMap(receivedFeedbackMap);     

    // master checks which ranks are connected to it
    for (auto &otherRank : receivedFeedbackMap) {
      for (auto &includedRank: otherRank.second) {
        if (utils::MasterSlave::_rank == includedRank) {
          _connectedRanks.push_back(otherRank.first);
          _mesh->getInitialCommunicationMap().push_back(otherRank.first);
        }
      }
    }    
  }
  else
  { // Slave
    utils::MasterSlave::_communication->broadcast(remoteParComSize, 0);
    for (int i=0; i < remoteParComSize; i++)
    {
      std::vector<int> initialFeedback;
      initialFeedback.push_back(-1);
      receivedFeedbackMap[i]=initialFeedback;
    }
    com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveFeedbackMap(receivedFeedbackMap);

    for (auto &otherRank : receivedFeedbackMap)
    {
      for (auto &includedRanks: otherRank.second)
      {
        if (utils::MasterSlave::_rank == includedRanks)
        {
          _connectedRanks.push_back(otherRank.first);
          _mesh->getInitialCommunicationMap().push_back(otherRank.first);
        }
      }
    }
  }
}

}}
