#pragma once

#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/Polation.hpp"

namespace precice {
namespace mapping {

/**
 * @brief Base class for interpolation based mappings, where mapping is done using a geometry-based linear combination of input values.
 *  Subclasses differ by the way computeMapping() fills the _interpolations and by mesh tagging. Mapping itself is shared.
 */
class BarycentricBaseMapping : public Mapping {
public:
  BarycentricBaseMapping(Constraint constraint, int dimensions);

  /// Destructor, empty.
  virtual ~BarycentricBaseMapping() = default;

  /// Computes the projections and interpolation relations. Must be done by inherited subclass
  virtual void computeMapping() = 0;

  virtual bool hasComputedMapping() const final;

  /// Removes a computed mapping.
  virtual void clear() final;

  /**
   * @brief Uses projection and interpolation relations to map data values.
   *
   * @pre computeMapping() has been called
   *
   * @param[in] inputDataID Data ID of input data values to be mapped from.
   * @param[in] outputDataID Data ID of output data values to be mapped to.
   */
  virtual void map(
      int inputDataID,
      int outputDataID) final;

  virtual void tagMeshFirstRound() final;
  virtual void tagMeshSecondRound() final;

private:
  logging::Logger _log{"mapping::BarycentricBaseMapping"};

protected:
  std::vector<Polation> _interpolations;

  bool _hasComputedMapping = false;
};

} // namespace mapping
} // namespace precice
