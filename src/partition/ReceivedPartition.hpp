#pragma once

#include <vector>
#include "Partition.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace partition {

/**
 * @brief A partition that is computed from a mesh received from another participant.
 *
 * A mesh is received by the master rank and re-partitioned among all slave ranks.
 * Afterwards necessary distribution data structures are set up.
 */
class ReceivedPartition : public Partition {
public:
  /// Defines the typ of geometric filter used
  enum GeometricFilter {
    /// undefined
    UNDEFINED,
    /// No geometric filter used (e.g. for RBF mappings)
    NO_FILTER,
    /// Filter at master and communicate only filtered mesh.
    FILTER_FIRST,
    /// Broadcast first and filter then
    BROADCAST_FILTER
  };

  /// Constructor
  ReceivedPartition(mesh::PtrMesh mesh, GeometricFilter geometricFilter, double safetyFactor);

  virtual ~ReceivedPartition() {}

  virtual void communicate() override;

  virtual void compute() override;

  virtual void communicateBoundingBox();
  virtual void computeBoundingBox();
  virtual void compareBoundingBoxes() override;

private:
  /// Sets _bb to the union with the mesh from fromMapping resp. toMapping, also enlage by _safetyFactor
  void prepareBoundingBox();

  /// Checks whether two bounding boxes are overlapping
  bool overlapping(mesh::Mesh::BoundingBox const &currentBB, mesh::Mesh::BoundingBox const &receivedBB);

  /// Checks if vertex in contained in _bb
  bool isVertexInBB(const mesh::Vertex &vertex);

  /** Checks whether provided meshes are empty.
   *
   * Empty provided meshes mean that the re-partitioning completely filtered
   * out the mesh received on this rank at the coupling interface.
   */
  bool areProvidedMeshesEmpty() const;

  virtual void createOwnerInformation() override;

  /// Helper function for 'createOwnerFunction' to set local owner information
  void setOwnerInformation(const std::vector<int> &ownerVec);

  GeometricFilter _geometricFilter;

  mesh::Mesh::BoundingBox _bb;

  int _dimensions;

  double _safetyFactor;

  logging::Logger _log{"partition::ReceivedPartition"};
};

} // namespace partition
} // namespace precice
