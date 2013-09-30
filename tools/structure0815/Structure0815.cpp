#include "Structure0815.hpp"
#include "utils/GeometryComputations.hpp"

Structure0815:: Structure0815
(
  int              dim,
  double           density,
  const DynVector& gravity,
  const DynVector& vertices,
  const tarch::la::DynamicVector<int> faces )
:
  _dim(dim),
  _density(density),
  _gravity(gravity),
  _centerOfGravity(dim, 0.0),
  _totalMass(0.0),
  _time(0.0),
  _timesteps(0),
  _vertices(vertices),
  _faces(faces),
  _forces(vertices.size(), 0.0),
  _velocities(vertices.size(), 0.0),
  _oldVelocities(vertices.size(), 0.0),
  _velocityDeltas(vertices.size(), 0.0),
  _displacements(vertices.size(), 0.0),
  _oldDisplacements(vertices.size(), 0.0),
  _displacementDeltas(vertices.size(), 0.0),
  _fixed(false),
  _fixture(dim, 0.0)
{
  double totalVolume = 0.0;
  DynVector zero(_dim, 0.0);
  if (_dim == 2){
    DynVector coords0(zero);
    DynVector coords1(zero);
    for (int iEdge=0; iEdge < _faces.size()/2; iEdge++){
      int index = iEdge * 2;
      for (int i=0; i<_dim; i++){
        coords0[i] = _vertices[_faces[index] * _dim + i];
        coords1[i] = _vertices[_faces[index+1] * _dim + i];
      }
      typedef precice::utils::GeometryComputations GeoComp;
      double area = GeoComp::triangleArea(zero, coords0, coords1);
      area = std::abs(area); // since it comes out signed from cross-prod
      totalVolume += area;
      if (not tarch::la::equals(area, 0.0)){
        _centerOfGravity += (coords0 + coords1) * area / 3.0;
      }
    }
  }
  else {
    assertion1(_dim == 3, _dim);
    DynVector coords0(zero);
    DynVector coords1(zero);
    DynVector coords2(zero);
    DynVector vec01(zero);
    DynVector vec02(zero);
    DynVector vec03(zero);
    DynVector crossVec(zero);
    for (int iTriangle=0; iTriangle < _faces.size()/3; iTriangle++){
      int index = iTriangle * 3;
      for (int i=0; i<_dim; i++){
        coords0[i] = _vertices[_faces[index] * _dim + i];
        coords1[i] = _vertices[_faces[index+1] * _dim + i];
        coords2[i] = _vertices[_faces[index+2] * _dim + i];
      }
      vec01 = coords1;
      vec01 -= coords0;
      vec02 = coords2;
      vec02 -= coords0;
      vec03 = zero;
      vec03 -= coords0;
      crossVec = tarch::la::cross(vec01, vec02, crossVec);
      //STRUCTURE_DEBUG("crossVector=" << crossVec);
      double volume = tarch::la::dot(crossVec, vec03) / 6.0;
      volume = std::abs(volume);
      //STRUCTURE_DEBUG("v0=" << coords0 << ", v1=" << coords1 << ", v2=" << coords2 << ", volume=" << volume);
      totalVolume += volume;
      if (not tarch::la::equals(volume, 0.0)){
        _centerOfGravity += (coords0 + coords1 + coords2) * volume / 4.0;
      }
    }
  }
  _totalMass = totalVolume * _density;
  _centerOfGravity /= totalVolume;
}

void Structure0815:: fixPoint
(
  const DynVector& fixture )
{
  _fixed = true;
  _fixture = fixture;
}

void Structure0815:: iterate
(
  double dt )
{
  DynVector zero(_dim, 0.0);
  DynVector force(3, 0.0);
  DynVector totalForce(_gravity);
  totalForce *= _totalMass; // Makes gravity acceleration a force
  DynVector r(3, 0.0); // Distance vector from center of gravity
  DynVector totalTorque(3, 0.0); // Always init with three components
  DynVector torque(3, 0.0);

  int verticesSize = _vertices.size() / _dim;
  for (int iVertex=0; iVertex < verticesSize; iVertex++){

    for (int i=0; i < _dim; i++){
      force[i] = _forces[iVertex*_dim + i];
      totalForce[i] += force[i];
    }

    for (int i=0; i<_dim; i++){
      r[i] = _vertices[iVertex*_dim + i] - _centerOfGravity[i];
    }
    tarch::la::cross(r, force, torque);
    totalTorque += torque;
  }

  STRUCTURE_DEBUG("Total mass = " << _totalMass);
  STRUCTURE_DEBUG("Total force = " << totalForce);
  STRUCTURE_DEBUG("Center of gravity = " << _centerOfGravity);
  STRUCTURE_DEBUG("Total torque = " << totalTorque);

  DynVector translVelocityDelta(zero);
  if (not _fixed){
    // Compute values of next timestep
    translVelocityDelta = (totalForce / _totalMass) * dt;
    STRUCTURE_DEBUG("translVelocityDelta = " << translVelocityDelta);
  }

  // Set values of next timestep
  DynVector rotForce(3, 0.0);
  DynVector rotVelocityDelta(zero);
  DynVector writeVelocity(zero);
  DynVector writeDisplacement(zero);
  DynVector velocityDelta(zero);
  double normR = 0.0;
  for (int iVertex=0; iVertex < verticesSize; iVertex++){
    int index = iVertex *_dim;
    // Compute and write velocity
    for (int i=0; i<_dim; i++){
      r[i] = _vertices[index + i] - _centerOfGravity[i];
    }
    normR = tarch::la::norm2(r);

    tarch::la::cross(totalTorque, r, rotForce);
    //rotForce[0] = torque * -1.0 * r[1]; // TODO
    //rotForce[1] = torque * r[0]; // TODO
    for (int i=0; i < _dim; i++){
      rotVelocityDelta[i] = (rotForce[i] / _totalMass) * dt;
    }
    velocityDelta = rotVelocityDelta + translVelocityDelta;
    for (int i=0; i < _dim; i++){
      _velocityDeltas[index+i] = velocityDelta[i];
      _velocities[index+i] = _oldVelocities[index+i] + _velocityDeltas[index+i];
      _displacementDeltas[index+i] =
          (_oldVelocities[index+i] + _velocityDeltas[index+i])*dt;
      _displacements[index+i] =
          _oldDisplacements[index+i] + _displacementDeltas[index+i];
    }
  }

  STRUCTURE_DEBUG("Computed time = " << _time
                  << ", computed timesteps = " << _timesteps);
}

void Structure0815:: timestep(double dt)
{
  STRUCTURE_DEBUG("Adding value changes to values.");
  _oldVelocities = _velocities;
  _oldDisplacements = _displacements;
  _time += dt;
  _timesteps++;
}

