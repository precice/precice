#pragma once

#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"

namespace precice {
namespace mapping {

/// Geometric multiscale mapping in radial direction
class RadialGeoMultiscaleMapping : public Mapping {
public:
  /**
   * @brief Geometric multiscale type of the mapping (spread or collect).
   *
   * A geometric multiscale mapping can either go from the 1D to the 2D/3D solver. Then, we call it "spread".
   * Or from the 2D/3D to the 1D solver, which we call "collect".
   */
  enum struct MultiscaleType {
    SPREAD,
    COLLECT
  };
  enum struct MultiscaleAxis {
    X = 0,
    Y = 1,
    Z = 2
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

  /// Takes care of compute-heavy operations needed only once to set up the mapping.
  void computeMapping() override;

  /// Removes a computed mapping.
  void clear() override;

  void tagMeshFirstRound() override;
  void tagMeshSecondRound() override;

  /// Returns name of the mapping
  std::string getName() const final override;

protected:
  /// @copydoc Mapping::mapConservative
  void mapConservative(const time::Sample &inData, Eigen::VectorXd &outData) override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData) override;

private:
  mutable logging::Logger _log{"mapping::RadialGeoMultiscaleMapping"};

  /// type of mapping, namely spread or collect
  MultiscaleType _type;

  /// main axis along which radial geometric multiscale coupling happens
  MultiscaleAxis _axis;

  /// computed vertex indices to map data from input vertices to output vertices and vice versa
  std::vector<int> _vertexIndicesSpread;
  std::vector<int> _vertexIndicesCollect;

  /// counts number of vertices between midpoints for averaging
  std::vector<int> _vertexCounter;
};

} // namespace mapping
} // namespace precice
