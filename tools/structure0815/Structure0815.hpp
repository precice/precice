#ifndef STRUCTURE0815_HPP_
#define STRUCTURE0815_HPP_

#include "Globals.hpp"
#include "io/TXTTableWriter.hpp"
#include "precice/SolverInterface.hpp"
#include "utils/Dimensions.hpp"

using Eigen::VectorXd;

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
class Structure0815 {
public:
  Structure0815(
      int                   dim,
      double                density,
      const VectorXd &      gravity,
      const VectorXd &      vertices,
      const Eigen::VectorXi faces);

  /// Only allows rotations around the point of fixture.
  void fixPoint(const VectorXd &fixture);

  /// Fixes translations in given directions.
  void fixTranslations(const Eigen::Matrix<bool, Eigen::Dynamic, 1> &fixedDirections);

  /**
   * @brief Sets the centerOfGravity and totalVolume manually.
   *
   * This can be useful, when the structure surface is shaped such that not
   * every node is visible from the origin.
   */
  void fixCharacteristics(
      const VectorXd &centerOfGravity,
      double          totalVolume);

  VectorXd &forces()
  {
    return _forces;
  }

  VectorXd &displacements()
  {
    return _displacements;
  }

  VectorXd &velocities()
  {
    return _velocities;
  }

  VectorXd &displacementDeltas()
  {
    return _displacementDeltas;
  }

  VectorXd &velocityDeltas()
  {
    return _velocityDeltas;
  }

  /// Performs an iteration from the last timestep.
  void iterate(double dt);

  /// Overwrites old timestep values with the current iteration.
  void timestep(double dt);

  const VectorXd &getCenterOfGravity() const
  {
    return _centerOfGravity;
  }

  double getMass() const
  {
    return _totalMass;
  }

private:
  // String constants for statistics writer entries
  const std::string TIMESTEPS;
  const std::string TIME;
  const std::string CENTEROFGRAVITY;
  const std::string TOTALMASS;
  const std::string TOTALVOLUME;

  int                                    _dim;
  double                                 _density;
  VectorXd                               _gravity;
  VectorXd                               _centerOfGravity;
  double                                 _totalMass;
  double                                 _time;
  int                                    _timesteps;
  VectorXd                               _vertices;
  Eigen::VectorXi                        _faces;
  VectorXd                               _forces;
  VectorXd                               _velocities;
  VectorXd                               _oldVelocities;
  VectorXd                               _velocityDeltas;
  VectorXd                               _displacements;
  VectorXd                               _oldDisplacements;
  VectorXd                               _displacementDeltas;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> _fixedTranslationDirections;
  bool                                   _fixed;
  VectorXd                               _fixture;
  bool                                   _fixedCharacteristics;
  precice::io::TXTTableWriter            _statisticsWriter;

  /**
   * @brief Computes the center of gravity, total mass and total volume.
   *
   * Does not store the values in the member variables.
   */
  void computeCharacteristics(
      VectorXd &centerOfGravity,
      double &  totalMass,
      double &  totalVolume);
};

#endif /* STRUCTURE0815_HPP_ */
