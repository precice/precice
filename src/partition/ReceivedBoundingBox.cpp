#include "partition/ReceivedBoundingBox.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/M2N.hpp"
#include "utils/EventTimings.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"
#include <vector>
#include <map>

using precice::utils::Event;

namespace precice {
namespace partition {

logging::Logger ReceivedBoundingBox:: _log ( "precice::partition::ReceivedBoundingBox" );

ReceivedBoundingBox::ReceivedBoundingBox
(
  mesh::PtrMesh mesh, double safetyFactor, GeometricFilter geometricFilter)
:
  Partition (mesh),
  _bb(mesh->getDimensions(), std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest())),
  _dimensions(mesh->getDimensions()),
  _safetyFactor(safetyFactor),
  _geometricFilter(geometricFilter)
{}

void ReceivedBoundingBox::communicateBoundingBox()
{  
  prepareBoundingBox();
  
  if (not utils::MasterSlave::_slaveMode)
  {
    _remoteParComSize=0;    
    _m2n->getMasterCommunication()->receive(_remoteParComSize, 0);

    // we need to construct _globalBB first and then we can fill in that
    mesh::Mesh::BoundingBox initialBB;
    for (int i=0; i < _dimensions; i++) {
      initialBB.push_back(std::make_pair(-1,-1));
    }    
    for (int remoteRank = 0; remoteRank < _remoteParComSize; remoteRank++ )
    {
      _globalBB[remoteRank]= initialBB;
    }

    // master receives global_bb from other master    
    com::CommunicateBoundingBox(_m2n->getMasterCommunication()).receiveBoundingBoxMap(_globalBB, 0 );
  }
}

void ReceivedBoundingBox::computeBoundingBox()
{
  // handle coupling mode first (i.e. serial participant)
  if (not utils::MasterSlave::_slaveMode && not utils::MasterSlave::_masterMode)
  { 
    //some code here
  }

  if (utils::MasterSlave::_slaveMode)
  {
    utils::MasterSlave::_communication->broadcast(_remoteParComSize, 0);

    // initializing the _globalBB
    for (int remoteRank = 0; remoteRank < _remoteParComSize; remoteRank++ )
    {
      _globalBB[remoteRank]= _bb;
    }
    
    // receive _globalBB from master
    com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveBoundingBoxMap(_globalBB);
    _numberOfRemoteRanks = _globalBB.size();      

    // received bounding boxes are comapred with each rank bb
    if (_numberOfRemoteRanks>0)
    {
      for (auto &other_rank: _globalBB)
      {
        if (compareBoundingBox(_bb,other_rank.second)) {
          _feedback.push_back(other_rank.first);            
        }
      }
    } 

    //send feedback to master
    utils::MasterSlave::_communication->send(_feedback, 0);     
  }
  else if (utils::MasterSlave::_masterMode)
  { // Master
    assertion(utils::MasterSlave::_rank==0);
    assertion(utils::MasterSlave::_size>1);
    // sending this participant communicator size to other master
    _m2n->getMasterCommunication()->send(utils::MasterSlave::_size , 0);
    utils::MasterSlave::_communication->broadcast(_remoteParComSize);
    com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendBoundingBoxMap(_globalBB);
    _numberOfRemoteRanks = _globalBB.size();

    // initializing _feedbackmap
    for (int rank_slave=1; rank_slave < utils::MasterSlave::_size ; rank_slave++)
    {
      _feedbackMap[rank_slave].push_back(-1);        
    }

    // master produces its own feedback
    for (auto &other_rank: _globalBB)
    {
      if (compareBoundingBox(_bb,other_rank.second))
      {
        _feedback.push_back(other_rank.first);
      }          
    }

    // master first adds ist feedback to feedbackmap
    _feedbackMap[0]=_feedback;

    // master receives feedbacks from slaves and put them into the feedbackmap
    for (int rank_slave=1; rank_slave < utils::MasterSlave::_size ; rank_slave++)
    {
      utils::MasterSlave::_communication->receive(_feedback, rank_slave);
      _feedbackMap[rank_slave]=_feedback;        
    }

    // master sends feedbackmap to other master
    com::CommunicateBoundingBox(_m2n->getMasterCommunication()).sendFeedbackMap(_feedbackMap,0);             
  }


  if (utils::MasterSlave::_masterMode)
  {
    for (int i = 0; i < _remoteParComSize; i++)
    {
      int counter=0;
      _m2n->getMasterCommunication()->receive(counter, 0);
      _vertexCounters.push_back(counter);
    }   
    utils::MasterSlave::_communication->broadcast(_vertexCounters);
  }
  else if (utils::MasterSlave::_slaveMode)
  {    
    utils::MasterSlave::_communication->broadcast(_vertexCounters, 0);  
  }  
}



}}
