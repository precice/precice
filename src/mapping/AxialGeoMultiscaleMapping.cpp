#include "AxialGeoMultiscaleMapping.hpp"
#include "mesh/Mesh.hpp"

namespace precice {

namespace mapping {

AxialGeoMultiscaleMapping::AxialGeoMultiscaleMapping(
    Constraint     constraint,
    int            dimensions,
    MultiscaleType type,
    AxialAxis      axis,
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
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
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
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
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

  int coord;
  if (_axis == X) {
    coord = 0;
  } else if (_axis == Y) {
    coord = 1;
  } else if (_axis == Z) {
    coord = 2;
  } else {
    PRECICE_ASSERT(false, "Unknown axis.");
  }

  PRECICE_ASSERT((inputValues.size() / inValueDimensions == static_cast<int>(input()->vertices().size())),
                 inputValues.size(), inValueDimensions, input()->vertices().size());
  PRECICE_ASSERT((outputValues.size() / outValueDimensions == static_cast<int>(output()->vertices().size())),
                 outputValues.size(), outValueDimensions, output()->vertices().size());

  PRECICE_DEBUG("Map consistent");
  if (_type == SPREAD) {
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
      PRECICE_ASSERT(((i * outValueDimensions) + coord) < outputValues.size(), ((i * outValueDimensions) + coord), outputValues.size())
      outputValues((i * outValueDimensions) + coord) = 2 * inputValues(0) * (1 - distance * distance);
    }
  } else {
    PRECICE_ASSERT(_type == COLLECT);
    PRECICE_ASSERT(output()->vertices().size() == 1);
    PRECICE_ASSERT(outputValues.size() == 1);
    outputValues(0)     = 0;
    size_t const inSize = input()->vertices().size();
    for (size_t i = 0; i < inSize; i++) {
      PRECICE_ASSERT(((i * inValueDimensions) + coord) < inputValues.size(), ((i * inValueDimensions) + coord), inputValues.size())
      outputValues(0) += inputValues((i * inValueDimensions) + coord);
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
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
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

} // namespace mapping
} // namespace precice
