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
  void onMappingComputed(mesh::PtrMesh origins, mesh::PtrMesh searchSpace) final override;

  /// name of the nng mapping
  std::string getName() const final override;

protected:
  /// @copydoc Mapping::mapConservative
  void mapConservative(const Sample &inData, Eigen::VectorXd &outData) final override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(const Sample &inData, Eigen::VectorXd &outData) final override;
};

} // namespace mapping
} // namespace precice
