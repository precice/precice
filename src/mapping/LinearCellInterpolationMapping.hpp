#pragma once

#include <list>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/BarycentricBaseMapping.hpp"
#include "mapping/Polation.hpp"

namespace precice {
namespace mapping {

/**
 * @brief Mapping using orthogonal projection to nearest triangle/edge/vertex and
 *        linear interpolation from projected point.
 */
class LinearCellInterpolationMapping : public BarycentricBaseMapping {
public:
  /// Constructor, taking mapping constraint.
  LinearCellInterpolationMapping(Constraint constraint, int dimensions);

  /// Computes the projections and interpolation relations.
  void computeMapping() override;

private:
  logging::Logger _log{"mapping::LinearCellInterpolationMappingMapping"};
};

} // namespace mapping
} // namespace precice
