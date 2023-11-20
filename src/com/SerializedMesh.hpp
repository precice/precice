#pragma once

#include <vector>

namespace precice {
namespace mesh {
class Mesh;
} // namespace mesh
namespace com {
class Communication;

namespace serialize {

/// serialized representation of mesh::Mesh
class SerializedMesh {
public:
  /** serializes a given mesh::Mesh
   *
   * Calls assertValid prior
   */
  static SerializedMesh serialize(const mesh::Mesh &mesh);

  /** adds the serialized mesh to an actual mesh
   *
   * The mesh to add to may contain vertices.
   */
  void addToMesh(mesh::Mesh &mesh) const;

  /// asserts the content for correctness
  void assertValid() const;

  void send(Communication &communication, int rankReceiver);

  /// receives a SerializedMesh and calls assertValid before returning
  static SerializedMesh receive(Communication &communication, int rankSender);

  void broadcastSend(Communication &communication);

  /// receives a SerializedMesh and calls assertValid before returning
  static SerializedMesh broadcastReceive(Communication &communication);

private:
  SerializedMesh() = default;

  /// contains the dimension, followed by the numbers of vertices, edges, triangles, and tetrahedra
  std::vector<int> sizes;

  /// sizes[0] * dimension coordinates for vertices
  std::vector<double> coords;

  // if no connectivity ( sum(sizes[1-3]) == 0)
  // then contains sizes[0] global IDs
  // else contains sizes[0] pairs of (global id, local id)
  //      followed by sizes[1] pairs of local ids defining edges
  //      followed by sizes[2] triples of local ids defining triangles
  //      followed by sizes[3] quadruples of local ids defining tetrahedra
  std::vector<int> ids;
};

} // namespace serialize
} // namespace com
} // namespace precice
