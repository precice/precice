#include "AxialGeoMultiscaleMapping.hpp"
#include "mesh/Mesh.hpp"

namespace precice::mapping {

AxialGeoMultiscaleMapping::AxialGeoMultiscaleMapping(
    Constraint     constraint,
    int            dimensions,
    MultiscaleType type,
    MultiscaleAxis axis,
    double         radius)
    : Mapping(constraint, dimensions),
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
    if (_type == SPREAD) {
      PRECICE_CHECK(input()->vertices().size() == 1, "You can only define an axial geometric multiscale mapping of type spread from a mesh with exactly one vertex.");
    } else {
      PRECICE_ASSERT(_type == COLLECT);
      PRECICE_CHECK(output()->vertices().size() == 1, "You can only define an axial geometric multiscale mapping of type spread from a mesh with exactly one vertex.");
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
  _hasComputedMapping = false;
}

void AxialGeoMultiscaleMapping::mapConservative(DataID inputDataID, DataID outputDataID)
{
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
  PRECICE_ASSERT(false, "Not yet implemented");
  PRECICE_DEBUG("Map conservative");
}

void AxialGeoMultiscaleMapping::mapConsistent(DataID inputDataID, DataID outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);

  const Eigen::VectorXd &inputValues        = input()->data(inputDataID)->values();
  int                    inValueDimensions  = input()->data(inputDataID)->getDimensions();
  int                    outValueDimensions = output()->data(outputDataID)->getDimensions();
  Eigen::VectorXd &      outputValues       = output()->data(outputDataID)->values();

  int effectiveCoordinate = _axis;
  PRECICE_ASSERT(effectiveCoordinate == 0 ||
                     effectiveCoordinate == 1 ||
                     effectiveCoordinate == 2,
                 "Unknown multiscale axis type.")

  PRECICE_ASSERT((inputValues.size() / inValueDimensions == static_cast<int>(input()->vertices().size())),
                 inputValues.size(), inValueDimensions, input()->vertices().size());
  PRECICE_ASSERT((outputValues.size() / outValueDimensions == static_cast<int>(output()->vertices().size())),
                 outputValues.size(), outValueDimensions, output()->vertices().size());

  PRECICE_DEBUG("Map consistent");
  if (_type == SPREAD) {
    /*
      3D vertices are assigned a value based on distance from the 1D vertex.
      Currently, a Hagen-Poiseuille profile determines the velocity value.
    */
    PRECICE_ASSERT(inputValues.size() == 1);
    PRECICE_ASSERT(input()->vertices().size() == 1);
    mesh::Vertex &v0      = input()->vertices()[0];
    size_t const  outSize = output()->vertices().size();

    for (size_t i = 0; i < outSize; i++) {
      Eigen::VectorXd difference(outValueDimensions);
      difference = v0.getCoords();
      difference -= output()->vertices()[i].getCoords();
      double distance = difference.norm() / _radius;
      PRECICE_CHECK(distance <= 1.05, "Output mesh has vertices that do not coincide with the geometric multiscale interface defined by the input mesh. Ratio of vertex distance to radius is {}.", distance);
      PRECICE_ASSERT(static_cast<int>((i * outValueDimensions) + effectiveCoordinate) < outputValues.size(), ((i * outValueDimensions) + effectiveCoordinate), outputValues.size())
      outputValues((i * outValueDimensions) + effectiveCoordinate) = 2 * inputValues(0) * (1 - distance * distance);
    }
  } else {
    PRECICE_ASSERT(_type == COLLECT);
    /*
      1D vertex is assigned the averaged value over all 3D vertices at the outlet,
      but only of the effectiveCoordinate component of the velocity vector.
    */
    PRECICE_ASSERT(output()->vertices().size() == 1);
    PRECICE_ASSERT(outputValues.size() == 1);
    outputValues(0)     = 0;
    size_t const inSize = input()->vertices().size();
    for (size_t i = 0; i < inSize; i++) {
      PRECICE_ASSERT(static_cast<int>((i * inValueDimensions) + effectiveCoordinate) < inputValues.size(), ((i * inValueDimensions) + effectiveCoordinate), inputValues.size())
      outputValues(0) += inputValues((i * inValueDimensions) + effectiveCoordinate);
    }
    outputValues(0) = outputValues(0) / inSize;
  }
}

void AxialGeoMultiscaleMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();

  computeMapping();

  if (getConstraint() == CONSISTENT) {
    PRECICE_ASSERT(_type == SPREAD, "Not yet implemented");
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
