#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/NearestNeighborBaseMapping.hpp"

namespace precice {
namespace mapping {

/// Mapping using nearest neighboring vertices and their local gradient values.
class NearestNeighborGradientMapping : public NearestNeighborBaseMapping {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   */
  NearestNeighborGradientMapping(Constraint constraint, int dimensions);

  /// Calculates the offsets needed for the gradient mappings after calculating the matched vertices
  void onMappingComputed(mesh::PtrMesh origins, mesh::PtrMesh searchSpace) override;

protected:
  /// @copydoc Mapping::mapConservative
  void mapConservative(DataID inputDataID, DataID outputDataID) override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(DataID inputDataID, DataID outputDataID) override;
};

} // namespace mapping
} // namespace precice
