#include "CommunicateBoundingBox.hpp"
#include "Communication.hpp"
#include "com/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include <map>
#include <vector>

namespace precice
{
namespace com
{

logging::Logger CommunicateBoundingBox::_log("precice::com::CommunicateBoundingBox");

CommunicateBoundingBox::CommunicateBoundingBox(
  com::PtrCommunication communication)
  : _communication(communication)
{}

void CommunicateBoundingBox::sendBoundingBox(
  const mesh::Mesh::BoundingBox &bb,
  int rankReceiver)
{
  TRACE(rankReceiver);

  for (auto &d : bb) {
    _communication->send(d.first, rankReceiver);
    _communication->send(d.second, rankReceiver);
  }
}

void CommunicateBoundingBox::receiveBoundingBox(
  mesh::Mesh::BoundingBox &bb,
  int rankSender)
{
  TRACE(rankSender);

  for (auto &d : bb) {
    _communication->receive(d.first, rankSender);
    _communication->receive(d.second, rankSender);
  }
}


void CommunicateBoundingBox::sendBoundingBoxMap(
  mesh::Mesh::BoundingBoxMap &bbm,
  int rankReceiver)
{
  TRACE(rankReceiver);

  for (auto &bb : bbm) {
      sendBoundingBox(bb.second, rankReceiver);
    }
  }


void CommunicateBoundingBox::receiveBoundingBoxMap(
  mesh::Mesh::BoundingBoxMap &bbm,
  int rankSender)
{
  TRACE(rankSender);

  for (auto &bb : bbm) {
    receiveBoundingBox(bb.second, rankSender);
  }
}

void CommunicateBoundingBox::broadcastSendBoundingBoxMap(
  mesh::Mesh::BoundingBoxMap &bbm)
{

  for (auto &rank : bbm) {
    for (auto &dimension : rank.second) {
      _communication->broadcast(dimension.first);
      _communication->broadcast(dimension.second);
    }
  }
  
}

void CommunicateBoundingBox::broadcastReceiveBoundingBoxMap(
  mesh::Mesh::BoundingBoxMap &bbm)
{
  for (auto &rank : bbm) {
    for (auto &dimension : rank.second) {
      _communication->broadcast(dimension.first, 0);
      _communication->broadcast(dimension.second, 0);
    }
  }
  
}


}
} // namespace precice, com
