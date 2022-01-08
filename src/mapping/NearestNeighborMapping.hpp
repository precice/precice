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
  NearestNeighborMapping(Constraint constraint, int dimensions);

  /// Maps input data to output data from input mesh to output mesh.
  virtual void map(int inputDataID, int outputDataID) override;

  /// Second lookup to calculate offsets. Only implemented for gradient mapping
  virtual void onMappingComputed(mesh::PtrMesh origins, mesh::PtrMesh searchSpace) override;
};

} // namespace mapping
} // namespace precice
