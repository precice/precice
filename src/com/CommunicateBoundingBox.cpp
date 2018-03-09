#include "CommunicateBoundingBox.hpp"
#include "Communication.hpp"
#include "com/SharedPointer.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
#include <map>
#include <vector>
#include <typeinfo>

namespace precice
{
namespace com
{

logging::Logger CommunicateBoundingBox::_log("precice::com::CommunicateBoundingBox");

CommunicateBoundingBox::CommunicateBoundingBox(
    com::PtrCommunication communication)
    : _communication(communication)
{
}

void CommunicateBoundingBox::sendBoundingBox(
    const mesh::Mesh::BoundingBox &bb,
    int                            rankReceiver)
{
  TRACE(rankReceiver);

  for (auto &d : bb) {
    _communication->send(d.first, rankReceiver);
    _communication->send(d.second, rankReceiver);
  }
}

void CommunicateBoundingBox::receiveBoundingBox(
    mesh::Mesh::BoundingBox &bb,
    int rankSender,
    int dim)
{
  TRACE(rankSender);
  
  assertion(typeid(bb)==typeid(mesh::Mesh::BoundingBox));
  assertion(bb.size()==dim);
  
  for (int d = 0; d < dim; d++) {
    _communication->receive(bb[d].first, rankSender);
    _communication->receive(bb[d].second, rankSender);
  }
}


void CommunicateBoundingBox::sendBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &gbb,
    int rankReceiver,
    int dim)
{
  TRACE(rankReceiver);

  int numOfRanks = gbb.size();
  _communication->send(numOfRanks, rankReceiver);

  for (auto &rank : gbb) {
    
    sendBoundingBox(rank.second, rankReceiver);
    
    }
  }


void CommunicateBoundingBox::receiveBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &gbb,
    int rankSender,
    int dim)
{
  TRACE(rankSender);

  int numOfRanks = gbb.size(); 
  mesh::Mesh::BoundingBox BoundingBox;

  _communication->receive(numOfRanks, rankSender);

  for (int i = 0; i < numOfRanks; i++) {    
    for (int d = 0; d < dim; d++){
      int min, max;
    _communication->receive(min, rankSender);
    _communication->receive(max, rankSender);
    BoundingBox.push_back(std::make_pair(min,max));
    }
    gbb[i] = BoundingBox;
    BoundingBox.clear();
  }
}


void CommunicateBoundingBox::receiveBoundingBoxFeedBack(
  std::map<int, std::vector<int>> &feedbackMap,
  int rankSender)
{

   //feedback members are sent to master
  int numberOfConnectedRanks;
  _communication->receive(numberOfConnectedRanks, rankSender);

  int connectedRank;
  
  for (int j=0; j<numberOfConnectedRanks; j++ )
  {
    _communication->receive(connectedRank, rankSender);
    feedbackMap[rankSender].push_back(connectedRank);
  }
}

void CommunicateBoundingBox::broadcastSendBoundingBox(
  mesh::Mesh::BoundingBoxMap &gbb,
  int dim)
{

  int numOfRanks = gbb.size();
  _communication->broadcast(numOfRanks);
  _communication->broadcast(dim);

  for (auto &rank : gbb) {
    for (auto &dimension : rank.second) {
    _communication->broadcast(dimension.first);
    _communication->broadcast(dimension.second);
    }
  }
}

void CommunicateBoundingBox::broadcastReceiveBoundingBox(
  mesh::Mesh::BoundingBoxMap &gbb,
  int dim)
{
  
  int numOfRanks;  
  double min, max;

  _communication->broadcast(&numOfRanks, 0);
  _communication->broadcast(&dim, 0);

  for (int rank=0; rank < numOfRanks; rank++) {
    for (int dimension=0; dimension<dim; dimension++) {
      _communication->broadcast(&min, 0);
      _communication->broadcast(&max, 0);
      gbb[rank].push_back(std::make_pair(min,max));
    }
  } 
}


}
} // namespace precice, com
