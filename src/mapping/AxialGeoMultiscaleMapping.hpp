#pragma once

#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"

namespace precice::mapping {

/// Geometric multiscale mapping in axial direction
class AxialGeoMultiscaleMapping : public Mapping {
public:
  /**
   * @brief Geometric multiscale nature of the mapping (spread or collect).
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
   * @brief Profile to use when type == SPREAD.
   * - UNIFORM: Uniform profile u(r) = U_mean
   * - PARABOLIC: laminar Poiseuille profile u(r) = 2*U_mean*(1 - (r/R)^2)
   */
  enum struct SpreadProfile {
    UNIFORM,
    PARABOLIC
  };
  enum struct MultiscaleDimension {
    D1D3,
    D1D2,
    D2D3
  };
  enum struct MultiscaleCrossSection {
    CIRCLE,
    SQUARE
  };
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   * @param[in] dimension Dimensionality pairing of the mapping
   * @param[in] type Geometric multiscale type of the mapping
   * @param[in] axis Main axis along which axial geometric multiscale coupling happens
   * @param[in] radius Radius of the 2D or 3D solver "tube"
   *            - For crossSection == CIRCLE: radius of the circular interface (R)
   *            - For crossSection == SQUARE: half side length of the square interface (side = 2*R)
   * @param[in] profile Profile for SPREAD (ignored for COLLECT).
   * @param[in] crossSection Shape of the interface cross-section (default: circle).
   */
  AxialGeoMultiscaleMapping(Constraint constraint, int dimensions, MultiscaleDimension dimension, MultiscaleType type, MultiscaleAxis axis, double radius, SpreadProfile profile = SpreadProfile::UNIFORM, MultiscaleCrossSection crossSection = MultiscaleCrossSection::CIRCLE);

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
  mutable logging::Logger _log{"mapping::AxialGeoMultiscaleMapping"};

  // dimensionality of mapping, namely 1D-3D, 1D-2D or 2D-3D
  MultiscaleDimension _dimension;

  /// type of mapping, namely spread or collect
  MultiscaleType _type;

  /// main axis along which axial geometric multiscale coupling happens
  MultiscaleAxis _axis;

  /// radius of the "tube" from or to which the data is mapped, i.e., radius of the circular interface between the two participants
  double _radius;

  /// selected profile used when _type == SPREAD
  SpreadProfile _profile;

  // cross-section of the pipe
  MultiscaleCrossSection _crossSection;

  /// computed vertex distances to map data from input vertex to output vertices
  std::vector<double> _vertexDistances;

  /// computed normalized transverse coordinates of oputput vertices for squared cross-section
  std::vector<Eigen::Vector2d> _vertexTransverseCoords;

  // Axis of 2D interface line
  int _lineCoord = -1;

  // nearest input vertex to output vertices
  std::vector<int> _nearestVertex;

  // band of nearest output vertices to input vertex
  std::vector<std::vector<int>> _collectBands;

  // Maximum distance of output vertices to each input vertex
  std::vector<double> _maxDistancePerInput;
};

} // namespace precice::mapping
