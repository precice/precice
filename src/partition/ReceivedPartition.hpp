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
 * A mesh is received by the primary rank and re-partitioned among all secondary ranks.
 * Afterwards necessary distribution data structures are set up.
 */
class ReceivedPartition : public Partition {
  /// Make the fixture friend of this class
  friend struct ReceivedPartitionFixture;

public:
  /// Defines the type of geometric filter used
  enum GeometricFilter {
    /// undefined
    UNDEFINED,
    /// No geometric filter used (e.g. for RBF mappings)
    NO_FILTER,
    /// Filter at primary rank and communicate only filtered mesh.
    ON_PRIMARY_RANK,
    /// Filter after communication on all secondary ranks
    ON_SECONDARY_RANKS
  };

  /// Constructor
  ReceivedPartition(const mesh::PtrMesh &mesh, GeometricFilter geometricFilter, double safetyFactor, bool allowDirectAccess = false);

  ~ReceivedPartition() override = default;

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

  /// Tag mesh in first round according to all mappings
  void tagMeshFirstRound();

  /// Tag mesh in second round according to all mappings
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

  bool _allowDirectAccess;

  logging::Logger _log{"partition::ReceivedPartition"};

  /// Max global vertex IDs of remote connected ranks
  std::vector<int> _remoteMaxGlobalVertexIDs;

  /// Min global vertex IDs of remote connected ranks
  std::vector<int> _remoteMinGlobalVertexIDs;
};

} // namespace partition
} // namespace precice
