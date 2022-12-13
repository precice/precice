#include <cstddef>
#include <memory>
#include <numeric>
#include <utility>

#include "CommunicateBoundingBox.hpp"
#include "Communication.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice::com::serialize {

class SerializedConnectionMap {
public:
  using ConnectionMap = std::map<Rank, std::vector<VertexID>>;

  static SerializedConnectionMap serialize(const ConnectionMap &cm)
  {
    SerializedConnectionMap scm;
    auto &                  content = scm.content;

    auto elements = std::transform_reduce(
        cm.begin(), cm.end(),
        1, // num entries
        std::plus<size_t>{}, [](const auto &kv) {
          // Rank, size, elements
          return 2 + kv.second.size();
        });
    content.reserve(elements);

    content.push_back(cm.size());
    for (const auto &[rank, ids] : cm) {
      content.push_back(rank);
      content.push_back(ids.size());
      content.insert(content.end(), ids.begin(), ids.end());
    }
    PRECICE_ASSERT(content.size() == static_cast<std::size_t>(elements));
    scm.assertValid();
    return scm;
  }

  ConnectionMap toConnectionMap() const
  {
    auto begin      = content.begin();
    int  numEntries = *begin;

    if (numEntries == 0) {
      return {};
    }
    std::advance(begin, 1);

    ConnectionMap cm;
    for (int entry = 0; entry < numEntries; ++entry) {
      // offset 0: rank
      auto rank = *begin;
      std::advance(begin, 1);
      // offset 1: vertexid count
      auto size = *begin;
      std::advance(begin, 1);
      // offset 2: n vertex ids
      auto verticesStart = begin;
      std::advance(begin, size);

      cm.emplace(rank, std::vector<int>(verticesStart, begin));
    }
    return cm;
  }

  void assertValid() const
  {
    const size_t totalSize = content.size();
    PRECICE_ASSERT(totalSize > 1);

    const auto numEntries = content.front();
    PRECICE_ASSERT(numEntries >= 0);

    ConnectionMap cm;
    size_t        consumed = 1;
    for (int entry = 0; entry < numEntries; ++entry) {
      // Minimal size is an empty entry: (rank, 0)
      PRECICE_ASSERT(totalSize >= consumed + 2);
      // offset 0: rank
      auto rank = content[consumed++];
      PRECICE_ASSERT(rank >= 0);
      // offset 1: vertexid count
      auto size = content[consumed++];
      PRECICE_ASSERT(size >= 0);

      // offset 2: n vertex ids
      PRECICE_ASSERT(totalSize >= consumed + size, "size larger than remaining");
      std::set<int> ids(content.begin() + consumed, content.begin() + consumed + size);
      PRECICE_ASSERT(ids.size() == static_cast<std::size_t>(size), "Duplicate vertex IDs");
      consumed += size;
    }
    PRECICE_ASSERT(consumed == totalSize, "Not everything consumed");
  }

  void sendTo(Communication &communication, int rankReceiver) const
  {
    communication.sendRange(content, rankReceiver);
  }

  static SerializedConnectionMap receiveFrom(Communication &communication, int rankSender)
  {
    SerializedConnectionMap scm;
    scm.content = communication.receiveRange(rankSender, AsVectorTag<int>{});
    scm.assertValid();
    return scm;
  }

  void broadcastSend(Communication &communication) const
  {
    communication.broadcast(content);
  }

  static SerializedConnectionMap broadcastReceive(Communication &communication)
  {
    SerializedConnectionMap scm;
    communication.broadcast(scm.content, 0);
    scm.assertValid();
    return scm;
  }

private:
  SerializedConnectionMap() = default;

  /// Num entries, Rank0, Size0, Entries0, Rank1, Size1, Entries0 ...
  /// @TODO move to size_t once we changed VertexIDs to size_t
  std::vector<int> content;
};

#if 0

class SerializedBoundingBoxMap {
public:
  using BoundingBoxMap = std::map<Rank, mesh::BoundingBox>;

  static SerializedBoundingBoxMap serialize(const BoundingBoxMap &bbm);
  BoundingBoxMap                  toBoundingBoxMap() const;

  void assertValid() const;

  void                            sendTo(Communication &communication, int rankReceiver);
  static SerializedBoundingBoxMap receiveFrom(Communication &communication, int rankSender);

  void                            broadcastSend(Communication &communication);
  static SerializedBoundingBoxMap broadcastReceive(Communication &communication);

private:
  SerializedBoundingBoxMap() = default;

  /// Num entries, Rank0, Rank1, ...
  std::vector<int> ranks;

  /** Dimensions followed by AABB coords
   * For 2D: 2, MinX0, MinY0, MaxX0, MaxY0, ...
   * For 3D: 3, MinX0, MinY0, MinZ0, MaxX0, MaxY0, MaxZ0, ...
   */
  std::vector<double> coords;
};

#endif
} // namespace precice::com::serialize

namespace precice::com {
CommunicateBoundingBox::CommunicateBoundingBox(
    com::PtrCommunication communication)
    : _communication(std::move(communication))
{
}

void CommunicateBoundingBox::sendBoundingBox(
    const mesh::BoundingBox &bb,
    int                      rankReceiver)
{
  PRECICE_TRACE(rankReceiver);
  _communication->sendRange(bb.dataVector(), rankReceiver);
}

void CommunicateBoundingBox::receiveBoundingBox(
    mesh::BoundingBox &bb,
    int                rankSender)
{
  PRECICE_TRACE(rankSender);
  auto              receivedData = _communication->receiveRange(rankSender, AsVectorTag<double>{});
  mesh::BoundingBox tempBB(receivedData);
  bb = std::move(tempBB);
}

void CommunicateBoundingBox::sendBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm,
    int                         rankReceiver)
{

  PRECICE_TRACE(rankReceiver);
  _communication->send(static_cast<int>(bbm.size()), rankReceiver);

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

  PRECICE_ASSERT(sizeOfReceivingMap == (int) bbm.size(), "Incoming size of map is not compatible");

  for (auto &bb : bbm) {
    receiveBoundingBox(bb.second, rankSender);
  }
}

void CommunicateBoundingBox::sendConnectionMap(
    std::map<int, std::vector<int>> const &fbm,
    int                                    rankReceiver)
{
  PRECICE_TRACE(rankReceiver);
  serialize::SerializedConnectionMap::serialize(fbm).sendTo(*_communication, rankReceiver);
}

///@todo needs some rewrite eventually. do we assume that the ranks are ordered or not? maybe change to vector
void CommunicateBoundingBox::receiveConnectionMap(
    std::map<int, std::vector<int>> &fbm,
    int                              rankSender)
{
  PRECICE_TRACE(rankSender);
  fbm = serialize::SerializedConnectionMap::receiveFrom(*_communication, rankSender).toConnectionMap();
}

void CommunicateBoundingBox::broadcastSendBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm)
{
  PRECICE_TRACE();
  _communication->broadcast(static_cast<int>(bbm.size()));

  for (const auto &rank : bbm) {
    _communication->broadcast(rank.second.dataVector());
  }
}

void CommunicateBoundingBox::broadcastReceiveBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm)
{
  PRECICE_TRACE();
  int sizeOfReceivingMap;
  _communication->broadcast(sizeOfReceivingMap, 0);
  PRECICE_ASSERT(sizeOfReceivingMap == (int) bbm.size());

  std::vector<double> receivedData;

  for (int i = 0; i < sizeOfReceivingMap; ++i) {
    _communication->broadcast(receivedData, 0);
    mesh::BoundingBox tempBB(receivedData);
    bbm.at(i) = std::move(tempBB);
  }
}

void CommunicateBoundingBox::broadcastSendConnectionMap(
    std::map<int, std::vector<int>> const &fbm)
{
  PRECICE_TRACE();
  serialize::SerializedConnectionMap::serialize(fbm).broadcastSend(*_communication);
}

void CommunicateBoundingBox::broadcastReceiveConnectionMap(
    std::map<int, std::vector<int>> &fbm)
{
  PRECICE_TRACE();
  fbm = serialize::SerializedConnectionMap::broadcastReceive(*_communication).toConnectionMap();
}

} // namespace precice::com
