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
class LinearCellInterpolation : public BarycentricBaseMapping {
public:
  /// Constructor, taking mapping constraint.
  LinearCellInterpolation(Constraint constraint, int dimensions);


  /// Computes the projections and interpolation relations.
  void computeMapping() override;

private:
  logging::Logger _log{"mapping::LinearCellInterpolationMapping"};
};

} // namespace mapping
} // namespace precice
