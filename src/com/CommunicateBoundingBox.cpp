#include "CommunicateBoundingBox.hpp"
#include "Communication.hpp"



namespace precice
{
namespace com
{
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

  for (const auto &d : bb) {
    _communication->send(d.first, rankReceiver);
    _communication->send(d.second, rankReceiver);
  }
}

void CommunicateBoundingBox::receiveBoundingBox(
    mesh::Mesh::BoundingBox &bb,
    int                      rankSender)
{
  TRACE(rankSender);

  for (auto &d : bb) {
    _communication->receive(d.first, rankSender);
    _communication->receive(d.second, rankSender);
  }
}

void CommunicateBoundingBox::sendBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm,
    int                         rankReceiver)
{
  TRACE(rankReceiver);

  for (const auto &bb : bbm) {
    sendBoundingBox(bb.second, rankReceiver);
  }
}

void CommunicateBoundingBox::receiveBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm,
    int                         rankSender)
{
  TRACE(rankSender);

  for (auto &bb : bbm) {
    receiveBoundingBox(bb.second, rankSender);
  }
}

void CommunicateBoundingBox::sendFeedbackMap(
  std::map<int, std::vector<int>> &fbm,
  int                         rankReceiver)
{
  TRACE(rankReceiver);

  for (const auto &vect : fbm) {
   _communication->send(vect.first, rankReceiver);
   _communication->send(vect.second, rankReceiver);
  }  
}

void CommunicateBoundingBox::receiveFeedbackMap(
  std::map<int, std::vector<int>> &fbm,
  int                         rankSender)
{
  TRACE(rankSender);

  for (const auto &vect : fbm) {
    int rank;
    std::vector<int> connected_ranks;
    _communication->receive(rank, rankSender);
    _communication->receive(connected_ranks, rankSender);
    fbm[rank]=connected_ranks;
  }  
}

void CommunicateBoundingBox::broadcastSendBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm)
{

  for (const auto &rank : bbm) {
    for (const auto &dimension : rank.second) {
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

void CommunicateBoundingBox::broadcastSendFeedbackMap(
  std::map<int, std::vector<int>> &fbm)
{
  for (auto &rank : fbm) {
    _communication->broadcast(rank.second);
  }
}

void CommunicateBoundingBox::broadcastReceiveFeedbackMap(
  std::map<int, std::vector<int>> &fbm)
{
  for (auto &rank : fbm) {
    _communication->broadcast(rank.second,0);
  }
}

} // namespace com
} // namespace precice
