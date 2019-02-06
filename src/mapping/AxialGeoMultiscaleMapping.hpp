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
   * @brief Specifies additional constraints for a mapping.
   *
   * A consistent mapping retains mean values. When mapping displacements, e.g.
   * rigid body motions are retained. A conservative mapping retains the sum of
   * the values. Values integrated over some area should be mapped conservative,
   * while area independent values such as pressure or stresses should be mapped
   * consistent.
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

  double _radius;
};

}} // namespace precice, mapping
