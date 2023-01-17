#pragma once

#include <map>
#include <vector>

#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"

namespace precice::com {
class Communication;

namespace serialize {

class SerializedConnectionMap {
public:
  using ConnectionMap = std::map<Rank, std::vector<VertexID>>;

  static SerializedConnectionMap serialize(const ConnectionMap &cm);

  ConnectionMap toConnectionMap() const;

  void assertValid() const;

  void send(Communication &communication, int rankReceiver) const;

  static SerializedConnectionMap receive(Communication &communication, int rankSender);

  void broadcastSend(Communication &communication) const;

  static SerializedConnectionMap broadcastReceive(Communication &communication);

private:
  SerializedConnectionMap() = default;

  /// Num entries, Rank0, Size0, Entries0, Rank1, Size1, Entries0 ...
  /// @TODO move to size_t once we changed VertexIDs to size_t
  std::vector<int> content;
};

class SerializedBoundingBox {
public:
  static SerializedBoundingBox serialize(const mesh::BoundingBox &bbm);

  mesh::BoundingBox toBoundingBox() const;

  void assertValid() const;

  void send(Communication &communication, int rankReceiver);

  static SerializedBoundingBox receive(Communication &communication, int rankSender);

private:
  SerializedBoundingBox() = default;

  /** AABB coords
   * For 2D: 2, MinX, MinY, MaxX, MaxY
   * For 3D: 3, MinX, MinY, MinZ, MaxX, MaxY, MaxZ
   */
  std::vector<double> coords;
};

class SerializedBoundingBoxMap {
public:
  using BoundingBoxMap = std::map<Rank, mesh::BoundingBox>;

  static SerializedBoundingBoxMap serialize(const BoundingBoxMap &bbm);

  BoundingBoxMap toBoundingBoxMap() const;

  void assertValid() const;

  void send(Communication &communication, int rankReceiver);

  static SerializedBoundingBoxMap receive(Communication &communication, int rankSender);

  void broadcastSend(Communication &communication);

  static SerializedBoundingBoxMap broadcastReceive(Communication &communication);

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

} // namespace serialize

} // namespace precice::com
