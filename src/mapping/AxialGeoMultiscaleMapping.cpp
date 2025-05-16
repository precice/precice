#include "AxialGeoMultiscaleMapping.hpp"
#include <Eigen/Core>
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
  PRECICE_TRACE(output()->nVertices());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
    if (_type == MultiscaleType::SPREAD) {
      PRECICE_CHECK(input()->nVertices() == 1, "You can only define an axial geometric multiscale mapping of type spread from a mesh with exactly one vertex.");

      /* When we add support for 1D meshes (https://github.com/precice/precice/issues/1669),
        we should check for valid dimensions combination.

        const int inDataDimensions  = input()->getDimensions();
        const int outDataDimensions = output()->getDimensions();

        PRECICE_CHECK(input()->getDimensions() == 1, "The input mesh on an axial geometric multiscale mapping can only be 1D at the moment, but it was defined to be {}.", input()->getDimensions())
        PRECICE_CHECK(output()->getDimensions() == 3, "The output mesh on an axial geometric multiscale mapping can only be 3D at the moment, but it was defined to be {}.", input()->getDimensions())
      */
      const int outDataDimensions = 3;

      int effectiveCoordinate = static_cast<std::underlying_type_t<MultiscaleType>>(_axis); // Convert enum struct to int
      PRECICE_ASSERT(effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::X) ||
                         effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Y) ||
                         effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Z),
                     "Unknown multiscale axis type.");

      // compute distances between 1D vertex and 3D vertices
      mesh::Vertex    &v0                           = input()->vertex(0);
      size_t const     outSize                      = output()->nVertices();
      constexpr double distance_to_radius_threshold = 1.05;

      _vertexDistances.clear();
      _vertexDistances.reserve(output()->nVertices());

      for (size_t i = 0; i < outSize; i++) {
        Eigen::VectorXd difference(outDataDimensions);
        difference = v0.getCoords();
        difference -= output()->vertex(i).getCoords();
        double distance_to_radius = difference.norm() / _radius;
        PRECICE_CHECK(distance_to_radius <= distance_to_radius_threshold, "Output mesh has vertices that do not coincide with the geometric multiscale interface defined by the input mesh. Ratio of vertex distance to radius is {} (which is larger than the assumed threshold of distance_to_radius_threshold).", distance_to_radius);
        _vertexDistances.push_back(distance_to_radius);
      }
    } else {
      PRECICE_ASSERT(_type == MultiscaleType::COLLECT);
      PRECICE_CHECK(output()->nVertices() == 1, "You can only define an axial geometric multiscale mapping of type spread from a mesh with exactly one vertex.");
      // Nothing to do here: A consistent collect mapping only averages all the values, independently of their locations, and this is done in the mapConsistent() method.
    }
  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_UNREACHABLE("Axial conservative geometric multiscale mapping is not implemented");
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
  PRECICE_UNREACHABLE("Axial conservative geometric multiscale mapping is not implemented");
  PRECICE_DEBUG("Map conservative");
}

void AxialGeoMultiscaleMapping::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();

  const int              inDataDimensions = inData.dataDims;
  const Eigen::VectorXd &inputValues      = inData.values;
  Eigen::VectorXd       &outputValues     = outData;
  // TODO: check if this needs to change when access to mesh dimension is possible
  const int outDataDimensions = outData.size() / output()->nVertices();

  int effectiveCoordinate = static_cast<std::underlying_type_t<MultiscaleType>>(_axis);
  PRECICE_ASSERT(effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::X) ||
                     effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Y) ||
                     effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Z),
                 "Unknown multiscale axis type.");

  // Check that the number of values for the input and output is right according to their dimensions
  PRECICE_ASSERT((inputValues.size() / static_cast<std::size_t>(inDataDimensions) == input()->nVertices()),
                 inputValues.size(), inDataDimensions, input()->nVertices());
  PRECICE_ASSERT((outputValues.size() / static_cast<std::size_t>(outDataDimensions) == output()->nVertices()),
                 outputValues.size(), outDataDimensions, output()->nVertices());

  // We currently don't support 1D data, so we need that the user specifies data of the same dimensions on both sides
  PRECICE_ASSERT(inDataDimensions == outDataDimensions);

  PRECICE_DEBUG("Map consistent");
  if (_type == MultiscaleType::SPREAD) {
    /*
      3D vertices are assigned a value based on distance from the 1D vertex.
      Currently, a Hagen-Poiseuille profile determines the velocity value.
    */
    PRECICE_ASSERT(input()->nVertices() == 1);
    size_t const outSize = output()->nVertices();

    for (size_t i = 0; i < outSize; i++) {
      PRECICE_ASSERT(static_cast<size_t>((i * outDataDimensions) + effectiveCoordinate) < static_cast<size_t>(outputValues.size()), ((i * outDataDimensions) + effectiveCoordinate), outputValues.size());
      // When adding support for 2D, remember that this should be 1.5 * inputValues(effectiveCoordinate) * (1 - (_vertexDistances[i] * _vertexDistances[i]));
      outputValues((i * outDataDimensions) + effectiveCoordinate) = 2 * inputValues(effectiveCoordinate) * (1 - (_vertexDistances[i] * _vertexDistances[i]));
    }
  } else {
    PRECICE_ASSERT(_type == MultiscaleType::COLLECT);
    /*
      1D vertex is assigned the averaged value over all 3D vertices at the outlet,
      but only of the effectiveCoordinate component of the velocity vector.
    */
    PRECICE_ASSERT(output()->nVertices() == 1);
    outputValues(effectiveCoordinate) = 0;
    size_t const inSize               = input()->nVertices();
    for (size_t i = 0; i < inSize; i++) {
      PRECICE_ASSERT(static_cast<size_t>((i * inDataDimensions) + effectiveCoordinate) < static_cast<size_t>(inputValues.size()), ((i * inDataDimensions) + effectiveCoordinate), inputValues.size());
      outputValues(effectiveCoordinate) += inputValues((i * inDataDimensions) + effectiveCoordinate);
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
    PRECICE_ASSERT(input()->nVertices() == 1);

    input()->tagAll();

  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_UNREACHABLE("Axial conservative geometric multiscale mapping is not implemented");
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
