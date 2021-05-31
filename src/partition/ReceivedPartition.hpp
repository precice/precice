#pragma once

#include <string>
#include <vector>
#include "Partition.hpp"
#include "logging/Logger.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace m2n {
class M2N;
} // namespace m2n

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
    ON_MASTER,
    /// Filter after communication on all slave ranks
    ON_SLAVES
  };

  /// Constructor
  ReceivedPartition(mesh::PtrMesh mesh, GeometricFilter geometricFilter, double safetyFactor);

  virtual ~ReceivedPartition() {}

  void communicate() override;

  void compute() override;

  void compareBoundingBoxes() override;

private:
  /// return the one m2n, a ReceivedPartition can only have one m2n
  m2n::M2N &m2n();

  void filterByBoundingBox();

  /// Sets _bb to the union with the mesh from fromMapping resp. toMapping, also enlage by _safetyFactor
  void prepareBoundingBox();

  /** Checks whether provided meshes are empty.
   *
   * Empty provided meshes mean that the re-partitioning completely filtered
   * out the mesh received on this rank at the coupling interface.
   */
  bool isAnyProvidedMeshNonEmpty() const;

  /// Returns whether any mapping is defined
  bool hasAnyMapping() const;

  /// Tag mesh in first round accoring to all mappings
  void tagMeshFirstRound();

  /// Tag mesh in second round accoring to all mappings
  void tagMeshSecondRound();

  void createOwnerInformation();

  /// Helper function for 'createOwnerFunction' to set local owner information
  void setOwnerInformation(const std::vector<int> &ownerVec);

  /// Is the local other (i.e. provided) bounding box already prepared (i.e. has prepareBoundingBox() been called)
  bool _boundingBoxPrepared = false;

  GeometricFilter _geometricFilter;

  mesh::BoundingBox _bb;

  int _dimensions;

  double _safetyFactor;

  logging::Logger _log{"partition::ReceivedPartition"};

  /// Max global vertex IDs of remote connected ranks
  std::vector<int> _remoteMaxGlobalVertexIDs;

  /// Min global vertex IDs of remote connected ranks
  std::vector<int> _remoteMinGlobalVertexIDs;
};

} // namespace partition
} // namespace precice
