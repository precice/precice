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
class NearestProjectionMapping : public BarycentricBaseMapping {
public:
  /// Constructor, taking mapping constraint.
  NearestProjectionMapping(Constraint constraint, int dimensions);

  /// Destructor, empty.
  virtual ~NearestProjectionMapping() {}

  /// Computes the projections and interpolation relations.
  virtual void computeMapping() override;


private:
  logging::Logger _log{"mapping::NearestNeighborProjectionMapping"};
};

} // namespace mapping
} // namespace precice
