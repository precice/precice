#include <cstddef>
#include <iterator>
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

class SerializedBoundingBoxMap {
public:
  using BoundingBoxMap = std::map<Rank, mesh::BoundingBox>;

  static SerializedBoundingBoxMap serialize(const BoundingBoxMap &bbm)
  {
    if (bbm.empty()) {
      SerializedBoundingBoxMap sbbm;
      sbbm.info = {0};
      return sbbm;
    }

    SerializedBoundingBoxMap sbbm;
    auto &                   info   = sbbm.info;
    auto &                   coords = sbbm.coords;

    // Entries, dimensions, ranks
    auto size = bbm.size();
    info.reserve(2 + size);
    const auto dims = bbm.begin()->second.getDimension();
    PRECICE_ASSERT(dims == 2 || dims == 3);
    info.push_back(size);
    info.push_back(dims);
    coords.reserve(size * dims);

    for (const auto &[rank, bb] : bbm) {
      info.push_back(rank);
      auto min = bb.minCorner();
      std::copy_n(min.data(), dims, std::back_inserter(coords));
      auto max = bb.maxCorner();
      std::copy_n(max.data(), dims, std::back_inserter(coords));
    }
    sbbm.assertValid();
    return sbbm;
  }

  BoundingBoxMap toBoundingBoxMap() const
  {
    if (info.size() == 1) {
      return {};
    }

    auto size = info[0];
    auto dims = info[1];

    BoundingBoxMap bbm;

    ///@todo replace the coord mess after refactoring of AABB to min and max points
    std::vector<double> buffer(dims * 2);

    auto rankIter  = std::next(info.begin(), 2);
    auto coordIter = coords.begin();
    for (int entry = 0; entry < size; ++entry) {
      // Copy coords into buffer
      // Input:  minX, minY, minZ, maxX, maxY, maxZ
      // Output: minX, maxX minY, maxY, minZ, maxZ
      for (int d = 0; d < dims; ++d) {
        auto offset        = d * 2;
        buffer[offset]     = coordIter[d];        // min
        buffer[offset + 1] = coordIter[d + dims]; // max
      }
      bbm.emplace(*rankIter, mesh::BoundingBox(buffer));

      // advance the input iterators
      std::advance(rankIter, 1);
      std::advance(coordIter, dims * 2);
    }
    return bbm;
  }

  void assertValid() const
  {
    PRECICE_ASSERT(!info.empty());
    if (info.size() == 1) {
      PRECICE_ASSERT(info.front() == 0);
      PRECICE_ASSERT(coords.empty());
      return;
    }

    auto numEntries = info[0];
    PRECICE_ASSERT(static_cast<size_t>(numEntries) + 2 == info.size());
    auto dims = info[1];
    PRECICE_ASSERT(dims == 2 || dims == 3);

    PRECICE_ASSERT(coords.size() == static_cast<size_t>(numEntries * dims * 2));
    for (int entry = 0; entry < numEntries; ++entry) {
      PRECICE_ASSERT(info[2 + entry] >= 0);
    }
  }

  void sendTo(Communication &communication, int rankReceiver)
  {
    communication.sendRange(info, rankReceiver);
    if (info.size() > 1) {
      communication.sendRange(coords, rankReceiver);
    }
  }

  static SerializedBoundingBoxMap receiveFrom(Communication &communication, int rankSender)
  {
    SerializedBoundingBoxMap sbbm;
    sbbm.info = communication.receiveRange(rankSender, AsVectorTag<int>{});
    if (sbbm.info.size() > 1) {
      sbbm.coords = communication.receiveRange(rankSender, AsVectorTag<double>{});
    }
    sbbm.assertValid();
    return sbbm;
  }

  void broadcastSend(Communication &communication)
  {
    communication.broadcast(info);
    if (info.size() > 1) {
      communication.broadcast(coords);
    }
  }

  static SerializedBoundingBoxMap broadcastReceive(Communication &communication)
  {
    SerializedBoundingBoxMap sbbm;
    communication.broadcast(sbbm.info, 0);
    if (sbbm.info.size() > 1) {
      communication.broadcast(sbbm.coords, 0);
    }
    sbbm.assertValid();
    return sbbm;
  }

private:
  SerializedBoundingBoxMap() = default;

  /** Num entries, Dimensions, Rank0, Rank1, ...
   *
   * If there are no entries, then the serialization cannot deduce the amount of dimensions.
   * In this case the serialization will only contain 0!
   */

  std::vector<int> info;

  /** AABB coords
   * For 2D: 2, MinX0, MinY0, MaxX0, MaxY0, ...
   * For 3D: 3, MinX0, MinY0, MinZ0, MaxX0, MaxY0, MaxZ0, ...
   */
  std::vector<double> coords;
};

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
  serialize::SerializedBoundingBoxMap::serialize(bbm).sendTo(*_communication, rankReceiver);
}

void CommunicateBoundingBox::receiveBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm,
    int                         rankSender)
{
  PRECICE_TRACE(rankSender);

  bbm = serialize::SerializedBoundingBoxMap::receiveFrom(*_communication, rankSender).toBoundingBoxMap();
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

  serialize::SerializedBoundingBoxMap::serialize(bbm).broadcastSend(*_communication);
}

void CommunicateBoundingBox::broadcastReceiveBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm)
{
  PRECICE_TRACE();

  bbm = serialize::SerializedBoundingBoxMap::broadcastReceive(*_communication).toBoundingBoxMap();
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
