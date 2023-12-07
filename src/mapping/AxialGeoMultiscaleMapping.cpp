#include "AxialGeoMultiscaleMapping.hpp"
#include <Eigen/src/Core/Matrix.h>
#include "mesh/Mesh.hpp"

namespace precice::mapping {

AxialGeoMultiscaleMapping::AxialGeoMultiscaleMapping(
    Constraint     constraint,
    int            dimensions,
    MultiscaleType type,
    MultiscaleAxis axis,
    double         radius)
    : Mapping(constraint, dimensions, false, Mapping::InitialGuessRequirement::None),
      _type(type),
      _axis(axis),
      _radius(radius)
{
  setInputRequirement(Mapping::MeshRequirement::VERTEX);
  setOutputRequirement(Mapping::MeshRequirement::VERTEX);
}

void AxialGeoMultiscaleMapping::computeMapping()
{
  PRECICE_TRACE(output()->vertices().size());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
    if (_type == MultiscaleType::SPREAD) {
      PRECICE_CHECK(input()->vertices().size() == 1, "You can only define an axial geometric multiscale mapping of type spread from a mesh with exactly one vertex.");

      /* When we add support for 1D meshes (https://github.com/precice/precice/issues/1669),
        we should check for valid dimensions combination.

        const int inValueDimensions  = input()->getDimensions();
        const int outValueDimensions = output()->getDimensions();

        PRECICE_CHECK(input()->getDimensions() == 1, "The input mesh on an axial geometric multiscale mapping can only be 1D at the moment, but it was defined to be {}.", input()->getDimensions())
        PRECICE_CHECK(output()->getDimensions() == 3, "The output mesh on an axial geometric multiscale mapping can only be 3D at the moment, but it was defined to be {}.", input()->getDimensions())
      */
      const int outValueDimensions = 3;

      int effectiveCoordinate = static_cast<std::underlying_type_t<MultiscaleType>>(_axis); // Convert enum struct to int
      PRECICE_ASSERT(effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::X) ||
                         effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Y) ||
                         effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Z),
                     "Unknown multiscale axis type.")

      // compute distances between 1D vertex and 3D vertices
      mesh::Vertex &v0      = input()->vertices()[0];
      size_t const  outSize = output()->vertices().size();

      for (size_t i = 0; i < outSize; i++) {
        Eigen::VectorXd difference(outValueDimensions);
        difference = v0.getCoords();
        difference -= output()->vertices()[i].getCoords();
        double distance = difference.norm() / _radius;
        PRECICE_CHECK(distance <= 1.05, "Output mesh has vertices that do not coincide with the geometric multiscale interface defined by the input mesh. Ratio of vertex distance to radius is {} (which is larger than the assumed threshold of 1.05).", distance);
        _vertexDistances.push_back(distance);
      }
    } else {
      PRECICE_ASSERT(_type == MultiscaleType::COLLECT);
      PRECICE_CHECK(output()->vertices().size() == 1, "You can only define an axial geometric multiscale mapping of type spread from a mesh with exactly one vertex.");
      // Nothing to do here: A consistent collect mapping only averages all the values, independently of their locations, and this is done in the mapConsistent() method.
    }
  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_ASSERT(false, "Not yet implemented");
    PRECICE_DEBUG("Compute conservative mapping");
  }
  _hasComputedMapping = true;
}

void AxialGeoMultiscaleMapping::clear()
{
  PRECICE_TRACE();
  _vertexDistances.clear();
  _hasComputedMapping = false;
}

void AxialGeoMultiscaleMapping::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
  PRECICE_ASSERT(false, "Not yet implemented");
  PRECICE_DEBUG("Map conservative");
}

void AxialGeoMultiscaleMapping::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();

  const int              inValueDimensions = inData.dataDims;
  const Eigen::VectorXd &inputValues       = inData.values;
  Eigen::VectorXd &      outputValues      = outData;
  // TODO: check if this needs to change when access to mesh dimension is possible
  const int outValueDimensions = outData.size() / output()->vertices().size();

  int effectiveCoordinate = static_cast<std::underlying_type_t<MultiscaleType>>(_axis);
  PRECICE_ASSERT(effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::X) ||
                     effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Y) ||
                     effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Z),
                 "Unknown multiscale axis type.")

  PRECICE_ASSERT((inputValues.size() / inValueDimensions == static_cast<int>(input()->vertices().size())),
                 inputValues.size(), inValueDimensions, input()->vertices().size());
  PRECICE_ASSERT((outputValues.size() / outValueDimensions == static_cast<int>(output()->vertices().size())),
                 outputValues.size(), outValueDimensions, output()->vertices().size());

  PRECICE_DEBUG("Map consistent");
  if (_type == MultiscaleType::SPREAD) {
    /*
      3D vertices are assigned a value based on distance from the 1D vertex.
      Currently, a Hagen-Poiseuille profile determines the velocity value.
    */
    PRECICE_ASSERT(input()->vertices().size() == 1);
    size_t const outSize = output()->vertices().size();

    for (size_t i = 0; i < outSize; i++) {
      PRECICE_ASSERT(static_cast<int>((i * outValueDimensions) + effectiveCoordinate) < outputValues.size(), ((i * outValueDimensions) + effectiveCoordinate), outputValues.size());
      // When adding support for 2D, remember that this should be 1.5 * inputValues(effectiveCoordinate) * (1 - (_vertexDistances[i] * _vertexDistances[i]));
      outputValues((i * outValueDimensions) + effectiveCoordinate) = 2 * inputValues(effectiveCoordinate) * (1 - (_vertexDistances[i] * _vertexDistances[i]));
    }
  } else {
    PRECICE_ASSERT(_type == MultiscaleType::COLLECT);
    /*
      1D vertex is assigned the averaged value over all 3D vertices at the outlet,
      but only of the effectiveCoordinate component of the velocity vector.
    */
    PRECICE_ASSERT(output()->vertices().size() == 1);
    outputValues(effectiveCoordinate) = 0;
    size_t const inSize               = input()->vertices().size();
    for (size_t i = 0; i < inSize; i++) {
      PRECICE_ASSERT(static_cast<int>((i * inValueDimensions) + effectiveCoordinate) < inputValues.size(), ((i * inValueDimensions) + effectiveCoordinate), inputValues.size())
      outputValues(effectiveCoordinate) += inputValues((i * inValueDimensions) + effectiveCoordinate);
    }
    outputValues(effectiveCoordinate) = outputValues(effectiveCoordinate) / inSize;
  }
}

void AxialGeoMultiscaleMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();

  computeMapping();

  if (getConstraint() == CONSISTENT) {
    PRECICE_ASSERT(_type == MultiscaleType::SPREAD, "Not yet implemented");
    PRECICE_ASSERT(input()->vertices().size() == 1);

    input()->vertices()[0].tag();

  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_ASSERT(false, "Not yet implemented");
  }

  clear();
}

void AxialGeoMultiscaleMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // no operation needed here for the moment
}

std::string AxialGeoMultiscaleMapping::getName() const
{
  return "axial-geomultiscale";
}

} // namespace precice::mapping
