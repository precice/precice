#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/NearestNeighborBaseMapping.hpp"

namespace precice {
namespace mapping {

/// Mapping using nearest neighboring vertices.
class NearestNeighborMapping : public NearestNeighborBaseMapping {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   */
  NearestNeighborMapping(Constraint constraint, int dimensions, Mapping::CouplingKind kind = Mapping::CouplingKind::SURFACE);

protected:
  /// @copydoc Mapping::mapConservative
  void mapConservative(DataID inputDataID, DataID outputDataID) override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(DataID inputDataID, DataID outputDataID) override;
};

} // namespace mapping
} // namespace precice
