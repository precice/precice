#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"

namespace precice {
namespace mapping {

/// Mapping using nearest neighboring vertices and (eventually) their local gradient values.
/// Base class for Nearest Neighbor Mapping and Nearest Neighbor Gradient
class NearestNeighborBaseMapping : public Mapping {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   */
  NearestNeighborBaseMapping(Constraint constraint, int dimensions, bool hasGradient, std::string mappingName,
                             std::string mappingNameShort);

  /// Destructor, empty.
  virtual ~NearestNeighborBaseMapping() = default;

  /// Checks if this is a Nearest Neighbor Gradient mapping.
  bool requireGradient();

  /// Computes the mapping coefficients from the in- and output mesh.
  void computeMapping() override;

  /// Returns true, if computeMapping() has been called.
  bool hasComputedMapping() const override;

  /// Removes a computed mapping.
  void clear() override;

  /** Matches the offsets needed for the gradient mapping
   * Does nothing by default
   */
  virtual void onMappingComputed(mesh::PtrMesh origins, mesh::PtrMesh searchSpace);

  void tagMeshFirstRound() override;
  void tagMeshSecondRound() override;

protected:
  /// NearestNeighborMapping or NearestNeighborGradientMapping
  std::string mappingName;

  /// nn or nng
  std::string mappingNameShort;

  mutable logging::Logger _log{"mapping::" + mappingName};

  /// Compute the vector offset between the matched vector and the source vector (needed for gradient mapping)
  /// Optimization: save this as an std::vector<double> and use an Eigen::Map to create an interface that uses the correct dimensions.
  std::vector<Eigen::VectorXd> _offsetsMatched;

  /// Computed output vertex indices to map data from input vertices to.
  std::vector<int> _vertexIndices;

private:
  /// Flag to indicate whether computeMapping() has been called.
  bool _hasComputedMapping = false;

  /// Flag if the mapping is a gradient mapping or not
  bool _requireGradient;
};

} // namespace mapping
} // namespace precice
