#ifndef STRUCTURE0815_HPP_
#define STRUCTURE0815_HPP_

#include "Globals.hpp"
#include "precice/SolverInterface.hpp"
#include "utils/Dimensions.hpp"
#include "io/TXTTableWriter.hpp"

using precice::utils::DynVector;

/**
 * @brief Rigid body motion structure solver for 2D and 3D.
 *
 * The solver geometry is derived from a surface polygon in 2D or surface
 * triangulation in 3D. The origin point (zero point) must be visible from every
 * surface node, in order to enable a correct computation of volume and mass.
 *
 * At the moment, no real rotations are used, but a simple implicit timestepping
 * scheme using the forces at every node is performed. Thus, the volume of the
 * structure can change during the simulation, despite being a rigid body motion
 * solver. A smaller timestep reduces the volume/mass error.
 *
 * Translations can be fixed in every direction, or a fixed point can be defined
 * around which (only) rotations can occur.
 *
 * Statistical data is written each timestep into a file
 * "structure0815-statistics.txt".
 */
class Structure0815
{
public:

  Structure0815(
    int              dim,
    double           density,
    const DynVector& gravity,
    const DynVector& vertices,
    const tarch::la::DynamicVector<int> faces );

  /**
   * @brief Only allows rotations around the point of fixture.
   */
  void fixPoint(const DynVector& fixture);

  /**
   * @brief Fixes translations in given directions.
   */
  void fixTranslations(const tarch::la::DynamicVector<bool>& fixedDirections);

  DynVector& forces() { return _forces; }

  DynVector& displacements() { return _displacements; }

  DynVector& velocities() { return _velocities; }

  DynVector& displacementDeltas() { return _displacementDeltas; }

  DynVector& velocityDeltas() { return _velocityDeltas; }

  /**
   * @brief Performs an iteration from the last timestep.
   */
  void iterate(double dt);

  /**
   * @brief Overwrites old timestep values with the current iteration.
   */
  void timestep(double dt);

  const DynVector& getCenterOfGravity() const { return _centerOfGravity; }

  double getMass() const { return _totalMass; }

private:

  // String constants for statistics writer entries
  const std::string TIMESTEPS;
  const std::string TIME;
  const std::string CENTEROFGRAVITY;
  const std::string TOTALMASS;
  const std::string TOTALVOLUME;

  int _dim;
  double _density;
  DynVector _gravity;
  DynVector _centerOfGravity;
  double _totalMass;
  double _time;
  int _timesteps;
  DynVector _vertices;
  tarch::la::DynamicVector<int> _faces;
  DynVector _forces;
  DynVector _velocities;
  DynVector _oldVelocities;
  DynVector _velocityDeltas;
  DynVector _displacements;
  DynVector _oldDisplacements;
  DynVector _displacementDeltas;
  tarch::la::DynamicVector<bool> _fixedTranslationDirections;
  bool _fixed;
  DynVector _fixture;
  precice::io::TXTTableWriter _statisticsWriter;

  /**
   * @brief Computes the center of gravity, total mass and total volume.
   *
   * Does not store the values in the member variables.
   */
  void computeCharacteristics (
    DynVector& centerOfGravity,
    double&    totalMass,
    double&    totalVolume );
};

#endif /* STRUCTURE0815_HPP_ */
