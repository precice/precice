#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"

namespace precice {
namespace mapping {

/// Mapping using nearest neighboring vertices.
class NearestNeighborMapping : public Mapping {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   */
  NearestNeighborMapping(Constraint constraint, int dimensions);

  /// Destructor, empty.
  virtual ~NearestNeighborMapping() {}

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() override;

  /// Returns true, if computeMapping() has been called.
  virtual bool hasComputedMapping() const override;

  /// Removes a computed mapping.
  virtual void clear() override;

  /// Maps input data to output data from input mesh to output mesh.
  virtual void map(
      int inputDataID,
      int outputDataID) override;

  virtual void tagMeshFirstRound() override;
  virtual void tagMeshSecondRound() override;

private:
  mutable logging::Logger _log{"mapping::NearestNeighborMapping"};

  /// Flag to indicate whether computeMapping() has been called.
  bool _hasComputedMapping = false;

  /// Computed output vertex indices to map data from input vertices to.
  std::vector<int> _vertexIndices;
};

} // namespace mapping
} // namespace precice
