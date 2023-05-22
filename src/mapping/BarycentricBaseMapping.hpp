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

  /// Removes a computed mapping.
  void clear() final override;

  void tagMeshFirstRound() final override;
  void tagMeshSecondRound() final override;

private:
  logging::Logger _log{"mapping::BarycentricBaseMapping"};

protected:
  /// @copydoc Mapping::mapConservative
  void mapConservative(const Sample &inData, Eigen::VectorXd &outData) override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(const Sample &inData, Eigen::VectorXd &outData) override;

  std::vector<Polation> _interpolations;
};

} // namespace mapping
} // namespace precice
