#pragma once

#include <string>
#include "Partition.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

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
  void communicate() override;

  /// All distribution data structures are set up.
  void compute() override;

  void compareBoundingBoxes() override;

private:
  void prepare();

  logging::Logger _log{"partition::ProvidedPartition"};
};

} // namespace partition
} // namespace precice
