#ifndef STRUCTURE0815_HPP_
#define STRUCTURE0815_HPP_

#include "Globals.hpp"
#include "precice/SolverInterface.hpp"
#include "utils/Dimensions.hpp"

using precice::utils::DynVector;

class Structure0815
{
public:

  Structure0815(
    int              dim,
    double           density,
    const DynVector& gravity,
    const DynVector& vertices,
    const tarch::la::DynamicVector<int> faces );

  void fixPoint(const DynVector& fixture);

  DynVector& forces() { return _forces; }

  DynVector& displacements() { return _displacements; }

  DynVector& velocities() { return _velocities; }

  DynVector& displacementDeltas() { return _displacementDeltas; }

  DynVector& velocityDeltas() { return _velocityDeltas; }

  void iterate(double dt);

  void timestep(double dt);

  const DynVector& getCenterOfGravity() const { return _centerOfGravity; }

  double getMass() const { return _totalMass; }

private:

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
  bool _fixed;
  DynVector _fixture;
};

#endif /* STRUCTURE0815_HPP_ */
