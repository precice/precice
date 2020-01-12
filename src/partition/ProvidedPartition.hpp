#pragma once

#include "Partition.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace partition {

/**
 * @brief A partition that is provided by the participant.
 *
 * The participant already provides a partition by calling setMeshVertices etc.
 * If required the mesh needs to be sent to another participant.
 * Furthermore, distribution data structures need to be set up.
 */
class ProvidedPartition : public Partition {
public:
  ProvidedPartition(mesh::PtrMesh mesh);

  virtual ~ProvidedPartition() {}

  /// The mesh is gathered and sent to another participant (if required)
  virtual void communicate() override;

  /// All distribution data structures are set up.
  virtual void compute() override;

  virtual void communicateBoundingBox();
  virtual void computeBoundingBox();

private:
  /// Sets owner=True on all vertices
  virtual void createOwnerInformation() override;

  logging::Logger _log{"partition::ProvidedPartition"};
};

} // namespace partition
} // namespace precice
