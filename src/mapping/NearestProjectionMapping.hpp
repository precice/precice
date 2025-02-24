#pragma once

#include "logging/Logger.hpp"
#include "mapping/BarycentricBaseMapping.hpp"
#include "mapping/Polation.hpp"

namespace precice::mapping {

/**
 * @brief Mapping using orthogonal projection to nearest triangle/edge/vertex and
 *        linear interpolation from projected point.
 */
class NearestProjectionMapping : public BarycentricBaseMapping {
public:
  /// Constructor, taking mapping constraint.
  NearestProjectionMapping(Constraint constraint, int dimensions);

  /// Computes the projections and interpolation relations.
  void computeMapping() final;

  /// name of the np mapping
  std::string getName() const final;

private:
  logging::Logger _log{"mapping::NearestNeighborProjectionMapping"};
};

} // namespace precice::mapping
