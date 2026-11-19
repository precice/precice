#pragma once

#include "MeshContext.hpp"
#include "partition/ProvidedPartition.hpp"

#include <memory>

namespace precice::impl {

/// Context for a mesh provided by this participant
struct ProvidedMeshContext : MeshContext {
  std::shared_ptr<precice::partition::ProvidedPartition> partition;
};

} // namespace precice::impl
