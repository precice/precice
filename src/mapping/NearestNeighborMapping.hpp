#pragma once

#include "mapping/Mapping.hpp"
#include "tarch/logging/Log.h"
#include <vector>

namespace precice {
namespace mapping {

/**
 * @brief Mapping using nearest neighboring vertices.
 */
class NearestNeighborMapping : public Mapping
{
public:

  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   */
  NearestNeighborMapping ( Constraint constraint, int dimensions );

  /// Destructor, empty.
  virtual ~NearestNeighborMapping() {}

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping();

  /// Returns true, if computeMapping() has been called.
  virtual bool hasComputedMapping() const;

  /// Removes a computed mapping.
  virtual void clear();

  /// Maps input data to output data from input mesh to output mesh.
  virtual void map (
    int inputDataID,
    int outputDataID );

  virtual bool doesVertexContribute(int vertexID) const;
  virtual bool isProjectionMapping() const;

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Flag to indicate whether computeMapping() has been called.
  bool _hasComputedMapping;

  // @brief Computed output vertex indices to map data from input vertices to.
  std::vector<int> _vertexIndices;
};

}} // namespace precice, mapping
