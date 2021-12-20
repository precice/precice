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
  virtual ~NearestNeighborBaseMapping() {}

  /// Checks if this is a Nearest Neighbor Gradient mapping.
  bool hasGradient();

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() override;

  /// Returns true, if computeMapping() has been called.
  virtual bool hasComputedMapping() const override;

  /// Removes a computed mapping.
  virtual void clear() override;

  /// Maps input data to output data from input mesh to output mesh.
  virtual void map(
      int inputDataID,
      int outputDataID) = 0;

  virtual void tagMeshFirstRound() override;
  virtual void tagMeshSecondRound() override;

protected:

  /// NearestNeighborMapping or NearestNeighborGradientMapping 
  std::string MAPPING_NAME;

  /// nn or nng
  std::string MAPPING_NAME_SHORT;

  mutable logging::Logger _log{"mapping::" + MAPPING_NAME};

  /// Compute the vector difference between the matched vector and the source vector (needed for gradient mapping)
  std::vector<Eigen::VectorXd> _distancesMatched;

  /// Computed output vertex indices to map data from input vertices to.
  std::vector<int> _vertexIndices;

private:
  /// Flag to indicate whether computeMapping() has been called.
  bool _hasComputedMapping = false;

  /// Flag if the mapping is a gradient mapping or not
  bool _hasGradient;



  
};

} // namespace mapping
} // namespace precice
