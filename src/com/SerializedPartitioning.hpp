#pragma once

#include <map>
#include <vector>

#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"

namespace precice::com {
class Communication;

namespace serialize {

/// serialized representation of ConnectionMap
class SerializedConnectionMap {
public:
  using ConnectionMap = std::map<Rank, std::vector<VertexID>>;

  /** serializes a given ConnectionMap
   *
   * Calls assertValid
   */
  static SerializedConnectionMap serialize(const ConnectionMap &cm);

  /// Builds and returns the connection map represented by the serialized state
  ConnectionMap toConnectionMap() const;

  /// asserts the content for correctness
  void assertValid() const;

  void send(Communication &communication, int rankReceiver) const;

  /// receives a SerializedConnectionMap and calls assertValid before returning
  static SerializedConnectionMap receive(Communication &communication, int rankSender);

  void broadcastSend(Communication &communication) const;

  /// receives a SerializedConnectionMap and calls assertValid before returning
  static SerializedConnectionMap broadcastReceive(Communication &communication);

private:
  SerializedConnectionMap() = default;

  /// Num entries, Rank0, Size0, Entries0, Rank1, Size1, Entries0 ...
  /// @TODO move to size_t once we changed VertexIDs to size_t
  std::vector<int> content;
};

/// serialized representation of a mesh::BoundingBox
class SerializedBoundingBox {
public:
  /** serializes a given mesh::BoundingBox
   *
   * Calls assertValid
   */
  static SerializedBoundingBox serialize(const mesh::BoundingBox &bbm);

  /// Builds and returns the BoundingBox represented by the serialized state
  mesh::BoundingBox toBoundingBox() const;

  /// asserts the content for correctness
  void assertValid() const;

  void send(Communication &communication, int rankReceiver);

  /// receives a SerializedBoundingBox and calls assertValid before returning
  static SerializedBoundingBox receive(Communication &communication, int rankSender);

private:
  SerializedBoundingBox() = default;

  /** AABB coords
   * For 2D: 2, MinX, MinY, MaxX, MaxY
   * For 3D: 3, MinX, MinY, MinZ, MaxX, MaxY, MaxZ
   */
  std::vector<double> coords;
};

/// serialized representation of a BoundingBoxMap
class SerializedBoundingBoxMap {
public:
  using BoundingBoxMap = std::map<Rank, mesh::BoundingBox>;

  /** serializes a given BoundingBoxMap
   *
   * Calls assertValid
   */
  static SerializedBoundingBoxMap serialize(const BoundingBoxMap &bbm);

  /// Builds and returns the BoundingBoxMap represented by the serialized state
  BoundingBoxMap toBoundingBoxMap() const;

  /// asserts the content for correctness
  void assertValid() const;

  void send(Communication &communication, int rankReceiver);

  /// receives a SerializedBoundingBoxMap and calls assertValid before returning
  static SerializedBoundingBoxMap receive(Communication &communication, int rankSender);

  void broadcastSend(Communication &communication);

  /// receives a SerializedBoundingBoxMap and calls assertValid before returning
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
