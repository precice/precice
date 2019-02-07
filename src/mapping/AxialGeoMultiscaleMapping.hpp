#pragma once

#include "mapping/Mapping.hpp"
#include "logging/Logger.hpp"
#include <vector>

namespace precice {
namespace mapping {

/// Geometric multiscale mapping in axial direction
class AxialGeoMultiscaleMapping : public Mapping
{
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

  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   * @param[in] type Geometric multiscale type of the mapping
   * @param[in] radius Radius of the 1D solver "tube"
   */
  AxialGeoMultiscaleMapping ( Constraint constraint, int dimensions, MultiscaleType type, double radius );

  /// Destructor, empty.
  virtual ~AxialGeoMultiscaleMapping() {}

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() override;

  /// Returns true, if computeMapping() has been called.
  virtual bool hasComputedMapping() const override;

  /// Removes a computed mapping.
  virtual void clear() override;

  /// Maps input data to output data from input mesh to output mesh.
  virtual void map (
    int inputDataID,
    int outputDataID ) override;

  virtual void tagMeshFirstRound() override;
  virtual void tagMeshSecondRound() override;

private:
  mutable logging::Logger _log{"mapping::AxialGeoMultiscaleMapping"};

  /// Flag to indicate whether computeMapping() has been called.
  bool _hasComputedMapping = false;

  MultiscaleType _type;

  /// radius of the 1D "tube" from or to which the data is mapped
  double _radius;

  /// scaling to make up for difference between max and avg value for a certain shape function
  double _scaling = 0.0;
};

}} // namespace precice, mapping
