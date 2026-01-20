#include <cstddef>
#include <iterator>
#include <memory>
#include <numeric>
#include <utility>

#include "com/Communication.hpp"
#include "com/SerializedPartitioning.hpp"
#include "utils/assertion.hpp"

namespace precice::com::serialize {

// SerializedConnectionMap

SerializedConnectionMap SerializedConnectionMap::serialize(const ConnectionMap &cm)
{
  SerializedConnectionMap scm;
  auto                   &content = scm.content;

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

SerializedConnectionMap::ConnectionMap SerializedConnectionMap::toConnectionMap() const
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

void SerializedConnectionMap::assertValid() const
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

void SerializedConnectionMap::send(Communication &communication, int rankReceiver) const
{
  communication.sendRange(content, rankReceiver);
}

SerializedConnectionMap SerializedConnectionMap::receive(Communication &communication, int rankSender)
{
  SerializedConnectionMap scm;
  scm.content = communication.receiveRange(rankSender, asVector<int>);
  scm.assertValid();
  return scm;
}

void SerializedConnectionMap::broadcastSend(Communication &communication) const
{
  communication.broadcast(content);
}

SerializedConnectionMap SerializedConnectionMap::broadcastReceive(Communication &communication)
{
  SerializedConnectionMap scm;
  communication.broadcast(scm.content, 0);
  scm.assertValid();
  return scm;
}

// SerializedBoundingBox

SerializedBoundingBox SerializedBoundingBox::serialize(const mesh::BoundingBox &bb)
{
  SerializedBoundingBox sbb;

  // Entries, dimensions, ranks
  const auto dims = bb.getDimension();
  PRECICE_ASSERT(dims == 2 || dims == 3);

  auto &coords = sbb.coords;
  coords.reserve(1 + 2 * dims);
  coords.push_back(dims);
  // copy point coordinates
  auto min = bb.minCorner();
  std::copy_n(min.data(), dims, std::back_inserter(coords));
  auto max = bb.maxCorner();
  std::copy_n(max.data(), dims, std::back_inserter(coords));

  sbb.assertValid();
  return sbb;
}

mesh::BoundingBox SerializedBoundingBox::toBoundingBox() const
{
  auto dims = static_cast<int>(coords.at(0));

  ///@todo replace the coord mess after refactoring of AABB to min and max points
  std::vector<double> buffer(dims * 2);

  // Copy coords into buffer
  // Input:  minX, minY, minZ, maxX, maxY, maxZ
  // Output: minX, maxX minY, maxY, minZ, maxZ
  for (int d = 0; d < dims; ++d) {
    auto offset        = d * 2;
    buffer[offset]     = coords[1 + d];        // min
    buffer[offset + 1] = coords[1 + d + dims]; // max
  }
  return mesh::BoundingBox(std::move(buffer));
}

void SerializedBoundingBox::assertValid() const
{
  PRECICE_ASSERT(!coords.empty());
  auto dims = static_cast<size_t>(coords.front());
  PRECICE_ASSERT(dims == 2 || dims == 3);
  PRECICE_ASSERT(coords.size() == 1 + dims * 2);
}

void SerializedBoundingBox::send(Communication &communication, int rankReceiver)
{
  communication.sendRange(coords, rankReceiver);
}

SerializedBoundingBox SerializedBoundingBox::receive(Communication &communication, int rankSender)
{
  SerializedBoundingBox sbb;
  sbb.coords = communication.receiveRange(rankSender, asVector<double>);
  sbb.assertValid();
  return sbb;
}

// SerializedBoundingBoxMap

SerializedBoundingBoxMap SerializedBoundingBoxMap::serialize(const BoundingBoxMap &bbm)
{
  if (bbm.empty()) {
    SerializedBoundingBoxMap sbbm;
    sbbm.info = {0};
    return sbbm;
  }

  SerializedBoundingBoxMap sbbm;
  auto                    &info   = sbbm.info;
  auto                    &coords = sbbm.coords;

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

SerializedBoundingBoxMap::BoundingBoxMap SerializedBoundingBoxMap::toBoundingBoxMap() const
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

void SerializedBoundingBoxMap::assertValid() const
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

void SerializedBoundingBoxMap::send(Communication &communication, int rankReceiver)
{
  communication.sendRange(info, rankReceiver);
  if (info.size() > 1) {
    communication.sendRange(coords, rankReceiver);
  }
}

SerializedBoundingBoxMap SerializedBoundingBoxMap::receive(Communication &communication, int rankSender)
{
  SerializedBoundingBoxMap sbbm;
  sbbm.info = communication.receiveRange(rankSender, asVector<int>);
  if (sbbm.info.size() > 1) {
    sbbm.coords = communication.receiveRange(rankSender, asVector<double>);
  }
  sbbm.assertValid();
  return sbbm;
}

void SerializedBoundingBoxMap::broadcastSend(Communication &communication)
{
  communication.broadcast(info);
  if (info.size() > 1) {
    communication.broadcast(coords);
  }
}

SerializedBoundingBoxMap SerializedBoundingBoxMap::broadcastReceive(Communication &communication)
{
  SerializedBoundingBoxMap sbbm;
  communication.broadcast(sbbm.info, 0);
  if (sbbm.info.size() > 1) {
    communication.broadcast(sbbm.coords, 0);
  }
  sbbm.assertValid();
  return sbbm;
}
} // namespace precice::com::serialize
