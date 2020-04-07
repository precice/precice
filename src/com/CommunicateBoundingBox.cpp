#include "CommunicateBoundingBox.hpp"
#include "Communication.hpp"
#include "mesh/BoundingBox.hpp"

namespace precice {
namespace com {
CommunicateBoundingBox::CommunicateBoundingBox(
    com::PtrCommunication communication)
    : _communication(communication)
{
}

void CommunicateBoundingBox::sendBoundingBox(
    const mesh::BoundingBox &bb,
    int                            rankReceiver)
{
  PRECICE_TRACE(rankReceiver);

  _communication->send(bb.data(), rankReceiver);

}

void CommunicateBoundingBox::receiveBoundingBox(
    mesh::BoundingBox &bb,
    int                      rankSender)
{
  std::vector<double> receivedData;
  PRECICE_TRACE(rankSender);
  _communication->receive(receivedData, rankSender);
  bb = mesh::BoundingBox::createFromData(receivedData);
}

void CommunicateBoundingBox::sendBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm,
    int                         rankReceiver)
{

  PRECICE_TRACE(rankReceiver);
  _communication->send((int) bbm.size(), rankReceiver);

  for (const auto &bb : bbm) {
    sendBoundingBox(bb.second, rankReceiver);
  }
}

void CommunicateBoundingBox::receiveBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm,
    int                         rankSender)
{
  PRECICE_TRACE(rankSender);
  int sizeOfReceivingMap;
  _communication->receive(sizeOfReceivingMap, rankSender);

  PRECICE_ASSERT(sizeOfReceivingMap == (int) bbm.size());

  for (auto &bb : bbm) {
    receiveBoundingBox(bb.second, rankSender);
  }
}

void CommunicateBoundingBox::sendConnectionMap(
    std::map<int, std::vector<int>> const &fbm,
    int                                    rankReceiver)
{
  PRECICE_TRACE(rankReceiver);
  _communication->send((int) fbm.size(), rankReceiver);

  for (const auto &vect : fbm) {
    _communication->send(vect.first, rankReceiver);
    _communication->send(vect.second, rankReceiver);
  }
}

///@todo needs some rewrite eventually. do we assume that the ranks are ordered or not? maybe change to vector
void CommunicateBoundingBox::receiveConnectionMap(
    std::map<int, std::vector<int>> &fbm,
    int                              rankSender)
{
  PRECICE_TRACE(rankSender);
  int sizeOfReceivingMap;
  _communication->receive(sizeOfReceivingMap, rankSender);
  PRECICE_ASSERT(sizeOfReceivingMap == (int) fbm.size());

  std::vector<int> connected_ranks;

  for (size_t i = 0; i < fbm.size(); ++i) {
    int rank;
    _communication->receive(rank, rankSender);
    _communication->receive(connected_ranks, rankSender);
    fbm[rank] = connected_ranks;
    connected_ranks.clear();
  }
}

void CommunicateBoundingBox::broadcastSendBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm)
{
  PRECICE_TRACE();
  _communication->broadcast((int) bbm.size());

  for (const auto &rank : bbm) {
    _communication->broadcast(rank.second.data());
  }
}

void CommunicateBoundingBox::broadcastReceiveBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm)
{
  PRECICE_TRACE();
  int sizeOfReceivingMap;
  _communication->broadcast(sizeOfReceivingMap, 0);
  PRECICE_ASSERT(sizeOfReceivingMap == (int) bbm.size());

  for (auto &rank : bbm) {
    std::vector<double> receivedData;
    _communication->receive(receivedData,0);
    rank.second = mesh::BoundingBox::createFromData(receivedData);
  }
}

void CommunicateBoundingBox::broadcastSendConnectionMap(
    std::map<int, std::vector<int>> const &fbm)
{
  PRECICE_TRACE();
  _communication->broadcast((int) fbm.size());

  for (auto &rank : fbm) {
    _communication->broadcast(rank.second);
  }
}

void CommunicateBoundingBox::broadcastReceiveConnectionMap(
    std::map<int, std::vector<int>> &fbm)
{
  PRECICE_TRACE();
  int sizeOfReceivingMap;
  _communication->broadcast(sizeOfReceivingMap, 0);
  PRECICE_ASSERT(sizeOfReceivingMap == (int) fbm.size());

  for (auto &rank : fbm) {
    _communication->broadcast(rank.second, 0);
  }
}

} // namespace com
} // namespace precice
