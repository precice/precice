#pragma once

#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"

namespace precice {
namespace mapping {

/// Geometric multiscale mapping in axial direction
class RadialGeoMultiscaleMapping : public Mapping {
public:
  /**
   * @brief Geometric multiscale nature of the mapping (spread or collect).
   *
   * A geometric multiscale mapping can either go from the 1D to the 2D/3D solver. Then, we call it "spread".
   * Or from the 2D/3D to the 1D solver, which we call "collect".
   */
  enum MultiscaleType {
    SPREAD,
    COLLECT
  };
  enum MultiscaleAxis {
    X,
    Y,
    Z
  };

  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   * @param[in] type Geometric multiscale type of the mapping
   * @param[in] axis Main axis along which radial geometric multiscale coupling happens
   */
  RadialGeoMultiscaleMapping(Constraint constraint, int dimensions, MultiscaleType type, MultiscaleAxis axis);

  /// Computes the mapping coefficients from the in- and output mesh.
  void computeMapping() override;

  /// Returns true, if computeMapping() has been called.
  // bool hasComputedMapping() const override; // check if needed at all

  /// Removes a computed mapping.
  void clear() override;

  void tagMeshFirstRound() override;
  void tagMeshSecondRound() override;

  /// Returns name of the mapping - TODO: needed for porting to develop
  // std::string getName() const final override;

protected:
  /// @copydoc Mapping::mapConservative
  void mapConservative(DataID inputDataID, DataID outputDataID) override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(DataID inputDataID, DataID outputDataID) override;

private:
  mutable logging::Logger _log{"mapping::RadialGeoMultiscaleMapping"};

  MultiscaleType _type;

  /// main axis along which radial geometric multiscale coupling happens
  MultiscaleAxis _axis;
};

} // namespace mapping
} // namespace precice