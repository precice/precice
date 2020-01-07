#include "Structure0815.hpp"

#include "math/geometry.hpp"
#include "math/math.hpp"

using Eigen::Vector3d;
using Eigen::VectorXd;

Structure0815::Structure0815(
    int                   dim,
    double                density,
    const VectorXd &      gravity,
    const VectorXd &      vertices,
    const Eigen::VectorXi faces)
    : TIMESTEPS("Timesteps"),
      TIME("Time"),
      CENTEROFGRAVITY("Center-of-gravity"),
      TOTALMASS("Total mass"),
      TOTALVOLUME("Total volume"),
      _dim(dim),
      _density(density),
      _gravity(gravity),
      _centerOfGravity(VectorXd::Constant(dim, 0.0)),
      _totalMass(0.0),
      _time(0.0),
      _timesteps(0),
      _vertices(vertices),
      _faces(faces),
      _forces(VectorXd::Zero(vertices.size())),
      _velocities(VectorXd::Zero(vertices.size())),
      _oldVelocities(VectorXd::Zero(vertices.size())),
      _velocityDeltas(VectorXd::Zero(vertices.size())),
      _displacements(VectorXd::Zero(vertices.size())),
      _oldDisplacements(VectorXd::Zero(vertices.size())),
      _displacementDeltas(VectorXd::Zero(vertices.size())),
      _fixedTranslationDirections(Eigen::Matrix<bool, Eigen::Dynamic, 1>::Constant(dim, false)),
      _fixed(false),
      _fixture(VectorXd::Constant(dim, 0.0)),
      _fixedCharacteristics(false),
      _statisticsWriter("structure0815-statistics.txt")
{
  double totalVolume = 0.0;
  computeCharacteristics(_centerOfGravity, _totalMass, totalVolume);
  STRUCTURE_INFO("Center of gravity: " << _centerOfGravity);
  STRUCTURE_INFO("Total mass: " << _totalMass);
  STRUCTURE_INFO("Total volume: " << totalVolume);

  _statisticsWriter.addData(TIMESTEPS, precice::io::TXTTableWriter::INT);
  _statisticsWriter.addData(TIME, precice::io::TXTTableWriter::DOUBLE);
  if (dim == 2) {
    _statisticsWriter.addData(CENTEROFGRAVITY, precice::io::TXTTableWriter::VECTOR2D);
  } else {
    _statisticsWriter.addData(CENTEROFGRAVITY, precice::io::TXTTableWriter::VECTOR3D);
  }
  _statisticsWriter.addData(TOTALMASS, precice::io::TXTTableWriter::DOUBLE);
  _statisticsWriter.addData(TOTALVOLUME, precice::io::TXTTableWriter::DOUBLE);

  _statisticsWriter.writeData(TIMESTEPS, _timesteps);
  _statisticsWriter.writeData(TIME, _time);
  if (dim == 2) {
    Eigen::Vector2d centerOfGravity(_centerOfGravity);
    _statisticsWriter.writeData(CENTEROFGRAVITY, centerOfGravity);
  } else {
    assertion(_dim == 3, _dim);
    Eigen::Vector3d centerOfGravity(_centerOfGravity);
    _statisticsWriter.writeData(CENTEROFGRAVITY, centerOfGravity);
  }
  _statisticsWriter.writeData(TOTALMASS, _totalMass);
  _statisticsWriter.writeData(TOTALVOLUME, totalVolume);
}

void Structure0815::fixPoint(
    const VectorXd &fixture)
{
  _fixed   = true;
  _fixture = fixture;
  STRUCTURE_INFO("Fixed point at: " << fixture << " (only rotations possible)");
}

void Structure0815::fixTranslations(
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> &fixedDirections)
{
  _fixedTranslationDirections = fixedDirections;
  STRUCTURE_INFO("Fixed translations: " << fixedDirections);
}

void Structure0815::fixCharacteristics(
    const VectorXd &centerOfGravity,
    double          totalVolume)
{
  _fixedCharacteristics = true;
  _centerOfGravity      = centerOfGravity;
  _totalMass            = totalVolume * _density;
  STRUCTURE_INFO("Fixed characteristics!");
  STRUCTURE_INFO("Center of gravity: " << _centerOfGravity);
  STRUCTURE_INFO("Total mass: " << _totalMass);
  STRUCTURE_INFO("Total volume: " << totalVolume);
}

void Structure0815::iterate(
    double dt)
{
  VectorXd zero  = VectorXd::Zero(_dim);
  Vector3d force = Vector3d::Zero();
  VectorXd totalForce(_gravity);
  totalForce *= _totalMass;                // Makes gravity acceleration a force
  Vector3d r           = Vector3d::Zero(); // Distance vector from center of gravity
  Vector3d totalTorque = Vector3d::Zero(); // Always init with three components
  Vector3d torque      = Vector3d::Zero();

  int verticesSize = _vertices.size() / _dim;
  for (int iVertex = 0; iVertex < verticesSize; iVertex++) {

    for (int i = 0; i < _dim; i++) {
      force[i] = _forces[iVertex * _dim + i];
      totalForce[i] += force[i];
    }

    for (int i = 0; i < _dim; i++) {
      r[i] = _vertices[iVertex * _dim + i] - _centerOfGravity[i];
    }
    torque = r.cross(force);
    totalTorque += torque;
  }

  //STRUCTURE_DEBUG("Total mass = " << _totalMass);
  STRUCTURE_DEBUG("Total force = " << totalForce);
  //STRUCTURE_DEBUG("Center of gravity = " << _centerOfGravity);
  STRUCTURE_DEBUG("Total torque = " << totalTorque);

  VectorXd translVelocityDelta(zero);
  if (not _fixed) {
    // Compute values of next timestep
    translVelocityDelta = (totalForce / _totalMass) * dt;
    STRUCTURE_DEBUG("translVelocityDelta = " << translVelocityDelta);
  }
  for (int i = 0; i < _dim; i++) {
    if (_fixedTranslationDirections[i]) {
      translVelocityDelta[i] = 0.0;
    }
  }

  // Set values of next timestep
  Vector3d rotForce = Vector3d::Zero();
  VectorXd rotVelocityDelta(zero);
  VectorXd writeVelocity(zero);
  VectorXd writeDisplacement(zero);
  VectorXd velocityDelta(zero);
  double   normR = 0.0;
  for (int iVertex = 0; iVertex < verticesSize; iVertex++) {
    int index = iVertex * _dim;
    // Compute and write velocity
    for (int i = 0; i < _dim; i++) {
      r[i] = _vertices[index + i] - _centerOfGravity[i];
    }
    normR = r.norm();

    rotForce = totalTorque.cross(r);
    //rotForce[0] = torque * -1.0 * r[1]; /// @todo
    //rotForce[1] = torque * r[0]; /// @todo
    for (int i = 0; i < _dim; i++) {
      rotVelocityDelta[i] = (rotForce[i] / _totalMass) * dt;
    }
    velocityDelta = rotVelocityDelta + translVelocityDelta;
    for (int i = 0; i < _dim; i++) {
      _velocityDeltas[index + i] = velocityDelta[i];
      _velocities[index + i]     = _oldVelocities[index + i] + _velocityDeltas[index + i];
      _displacementDeltas[index + i] =
          (_oldVelocities[index + i] + _velocityDeltas[index + i]) * dt;
      _displacements[index + i] =
          _oldDisplacements[index + i] + _displacementDeltas[index + i];
    }
  }

  STRUCTURE_DEBUG("Computed time = " << _time
                                     << ", computed timesteps = " << _timesteps);
}

void Structure0815::timestep(double dt)
{
  STRUCTURE_DEBUG("Adding value changes to values.");
  _oldVelocities    = _velocities;
  _oldDisplacements = _displacements;
  _time += dt;
  _timesteps++;

  VectorXd centerOfGravity = VectorXd::Zero(_dim);
  double   totalMass       = 0.0;
  double   totalVolume     = 0.0;
  if (not _fixedCharacteristics) {
    computeCharacteristics(centerOfGravity, totalMass, totalVolume);
    STRUCTURE_INFO("Center of gravity delta: " << _centerOfGravity - centerOfGravity);
    STRUCTURE_INFO("Total mass delta: " << _totalMass - totalMass);
    STRUCTURE_INFO("Total volume delta: " << _totalMass / _density - totalVolume);
  }

  _statisticsWriter.writeData(TIMESTEPS, _timesteps);
  _statisticsWriter.writeData(TIME, _time);
  if (_dim == 2) {
    Eigen::Vector2d centerOfGravity2D(centerOfGravity);
    _statisticsWriter.writeData(CENTEROFGRAVITY, centerOfGravity2D);
  } else {
    assertion(_dim == 3, _dim);
    Eigen::Vector3d centerOfGravity3D(centerOfGravity);
    _statisticsWriter.writeData(CENTEROFGRAVITY, centerOfGravity3D);
  }
  _statisticsWriter.writeData(TOTALMASS, totalMass);
  _statisticsWriter.writeData(TOTALVOLUME, totalVolume);
}

void Structure0815::computeCharacteristics(
    VectorXd &centerOfGravity,
    double &  totalMass,
    double &  totalVolume)
{
  VectorXd zero   = VectorXd::Zero(_dim);
  centerOfGravity = zero;
  totalMass       = 0.0;
  totalVolume     = 0.0;

  if (_dim == 2) {
    VectorXd coords0(zero);
    VectorXd coords1(zero);
    for (int iEdge = 0; iEdge < _faces.size() / 2; iEdge++) {
      int index = iEdge * 2;
      for (int i = 0; i < _dim; i++) {
        coords0[i] = _vertices[_faces[index] * _dim + i];
        coords0[i] += _displacements[_faces[index] * _dim + i];
        coords1[i] = _vertices[_faces[index + 1] * _dim + i];
        coords1[i] += _displacements[_faces[index + 1] * _dim + i];
      }
      double area = geometry::triangleArea(zero, coords0, coords1);
      area        = std::abs(area); // since it comes out signed from cross-prod
      totalVolume += area;
      if (not precice::math::equals(area, 0.0)) {
        centerOfGravity += (coords0 + coords1) * area / 3.0;
      }
    }
  } else {
    assertion(_dim == 3, _dim);
    VectorXd coords0(zero);
    VectorXd coords1(zero);
    VectorXd coords2(zero);
    Vector3d vec01(zero);
    Vector3d vec02(zero);
    VectorXd vec03(zero);
    Vector3d crossVec(zero);
    for (int iTriangle = 0; iTriangle < _faces.size() / 3; iTriangle++) {
      int index = iTriangle * 3;
      for (int i = 0; i < _dim; i++) {
        coords0[i] = _vertices[_faces[index] * _dim + i];
        coords0[i] += _displacements[_faces[index] * _dim + i];
        coords1[i] = _vertices[_faces[index + 1] * _dim + i];
        coords1[i] += _displacements[_faces[index + 1] * _dim + i];
        coords2[i] = _vertices[_faces[index + 2] * _dim + i];
        coords2[i] += _displacements[_faces[index + 2] * _dim + i];
      }
      vec01 = coords1;
      vec01 -= coords0;
      vec02 = coords2;
      vec02 -= coords0;
      vec03 = zero;
      vec03 -= coords0;
      crossVec      = vec01.cross(vec02);
      double volume = crossVec.dot(vec03) / 6.0;
      volume        = std::abs(volume);
      totalVolume += volume;
      if (not precice::math::equals(volume, 0.0)) {
        centerOfGravity += (coords0 + coords1 + coords2) * volume / 4.0;
      }
    }
  }
  totalMass = totalVolume * _density;
  centerOfGravity /= totalVolume;
}
