#pragma once

#include <list>
#include <string>
#include <vector>
#include "Mapping.hpp"
#include "logging/Logger.hpp"
#include "query/FindClosest.hpp"

namespace precice {
namespace mapping {

/**
 * @brief Mapping using orthogonal projection to nearest triangle/edge/vertex and
 *        linear interpolation from projected point.
 */
class NearestProjectionMapping : public Mapping {
public:
  /// Constructor, taking mapping constraint.
  NearestProjectionMapping(Constraint constraint, int dimensions);

  /// Destructor, empty.
  virtual ~NearestProjectionMapping() {}

  /// Computes the projections and interpolation relations.
  virtual void computeMapping() override;

  virtual bool hasComputedMapping() const override;

  /// Removes a computed mapping.
  virtual void clear() override;

  /**
   * @brief Uses projection and interpolation relations to map data values.
   *
   * Preconditions:
   * - computeMapping() has been called
   *
   * @param[in] inputDataID Data ID of input data values to be mapped from.
   * @param[in] outputDataID Data ID of output data values to be mapped to.
   */
  virtual void map(
      int inputDataID,
      int outputDataID) override;

  virtual void tagMeshFirstRound() override;
  virtual void tagMeshSecondRound() override;

private:
  logging::Logger _log{"mapping::NearestProjectionMapping"};

  using InterpolationElements = std::vector<query::InterpolationElement>;
  std::vector<InterpolationElements> _weights;

  bool _hasComputedMapping = false;
};

} // namespace mapping
} // namespace precice
