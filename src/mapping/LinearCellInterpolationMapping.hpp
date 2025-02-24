#pragma once

#include "logging/Logger.hpp"
#include "mapping/BarycentricBaseMapping.hpp"

namespace precice::mapping {

/**
 * @brief Mapping using orthogonal projection to nearest triangle/edge/vertex and
 *        linear interpolation from projected point.
 */
class LinearCellInterpolationMapping : public BarycentricBaseMapping {
public:
  /// Constructor, taking mapping constraint.
  LinearCellInterpolationMapping(Constraint constraint, int dimensions);

  /// Computes the projections and interpolation relations.
  void computeMapping() final;

  /// name of the lci mapping
  std::string getName() const final;

private:
  logging::Logger _log{"mapping::LinearCellInterpolationMappingMapping"};
};

} // namespace precice::mapping
