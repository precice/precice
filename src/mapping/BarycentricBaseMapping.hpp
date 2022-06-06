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

  /// Removes a computed mapping.
  virtual void clear() final;

  virtual void tagMeshFirstRound() final;
  virtual void tagMeshSecondRound() final;

private:
  logging::Logger _log{"mapping::BarycentricBaseMapping"};

protected:
  /// @copydoc Mapping::mapConservative
  void mapConservative(DataID inputDataID, DataID outputDataID) override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(DataID inputDataID, DataID outputDataID) override;

  std::vector<Polation> _interpolations;
};

} // namespace mapping
} // namespace precice
