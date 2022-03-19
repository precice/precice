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
class VolumeCellInterpolation : public BarycentricBaseMapping {
public:
  /// Constructor, taking mapping constraint.
  VolumeCellInterpolation(Constraint constraint, int dimensions);

  /// Destructor, empty.
  virtual ~VolumeCellInterpolation() {}

  /// Computes the projections and interpolation relations.
  virtual void computeMapping() override;

  virtual void tagMeshFirstRound() override;
  virtual void tagMeshSecondRound() override;
};

} // namespace mapping
} // namespace precice
