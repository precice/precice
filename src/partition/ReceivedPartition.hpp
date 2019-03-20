#pragma once

#include <vector>
#include "Partition.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"

namespace precice
{
namespace partition
{

/**
 * @brief A partition that is computed from a mesh received from another participant.
 *
 * A mesh is received by the master rank and re-partitioned among all slave ranks.
 * Afterwards necessary distribution data structures are set up.
 */
class ReceivedPartition : public Partition
{
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

private:
  /// Create filteredMesh from the filtered _mesh.
  /*
   * Copies all vertices/edges/triangles that are either contained in the bounding box
   * or tagged to the filteredMesh. Edges and triangles are copied, when ALL vertices
   * are part of the filteredMesh i.e. their IDs are contained in vertexMap.
   */
  void filterMesh(mesh::Mesh &filteredMesh, const bool filterByBB);
  
  /// Sets _bb to the union with the mesh from fromMapping resp. toMapping, also enlage by _safetyFactor
  void prepareBoundingBox();

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
