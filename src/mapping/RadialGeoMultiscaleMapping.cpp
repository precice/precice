#include "RadialGeoMultiscaleMapping.hpp"
#include "mesh/Mesh.hpp"

namespace precice::mapping {

RadialGeoMultiscaleMapping::RadialGeoMultiscaleMapping(
    Constraint     constraint,
    int            dimensions,
    MultiscaleType type,
    MultiscaleAxis axis)
    : Mapping(constraint, dimensions),
      _type(type), _axis(axis)
{
  setInputRequirement(Mapping::MeshRequirement::VERTEX);
  setOutputRequirement(Mapping::MeshRequirement::VERTEX);
}

void RadialGeoMultiscaleMapping::computeMapping()
{
  PRECICE_TRACE(output()->vertices().size());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_ASSERT(false, "Not yet implemented");
    PRECICE_DEBUG("Compute conservative mapping");
  }
  _hasComputedMapping = true;
}

void RadialGeoMultiscaleMapping::clear()
{
  PRECICE_TRACE();
  _hasComputedMapping = false;
}

void RadialGeoMultiscaleMapping::mapConservative(DataID inputDataID, DataID outputDataID)
{
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
  PRECICE_ASSERT(false, "Not yet implemented");
  PRECICE_DEBUG("Map conservative");
}

void RadialGeoMultiscaleMapping::mapConsistent(DataID inputDataID, DataID outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);

  const Eigen::VectorXd &inputValues     = input()->data(inputDataID)->values();
  int                    valueDimensions = input()->data(inputDataID)->getDimensions();
  Eigen::VectorXd &      outputValues    = output()->data(outputDataID)->values();

  int effectiveCoordinate = _axis;
  PRECICE_ASSERT(effectiveCoordinate == 0 ||
                     effectiveCoordinate == 1 ||
                     effectiveCoordinate == 2,
                 "Unknown multiscale axis type.")

  PRECICE_ASSERT((inputValues.size() / valueDimensions == static_cast<int>(input()->vertices().size())),
                 inputValues.size(), valueDimensions, input()->vertices().size());
  PRECICE_ASSERT((outputValues.size() / valueDimensions == static_cast<int>(output()->vertices().size())),
                 outputValues.size(), valueDimensions, output()->vertices().size());
  PRECICE_ASSERT(valueDimensions == 1);

  PRECICE_DEBUG("Map consistent");
  if (_type == SPREAD) {
    /*
      3D vertices are projected onto the 1D axis and the data is then mapped
      to the nearest neighbors of the 1D vertices in projection space.
    */
    size_t const inSize  = input()->vertices().size();
    size_t const outSize = output()->vertices().size();

    // determine principle axis midpoints as borders to assign 3D vertices to their respective 1D vertex
    Eigen::VectorXd axisMidpoints(inSize);
    for (size_t i = 0; i < (inSize - 1); i++) {
      auto axisPositionCurrent = input()->vertices()[i].getCoords()[effectiveCoordinate];
      auto axisPositionNext    = input()->vertices()[i + 1].getCoords()[effectiveCoordinate];
      axisMidpoints(i)         = (axisPositionCurrent + axisPositionNext) / 2;
    }
    axisMidpoints(inSize - 1) = std::numeric_limits<double>::max(); // large number, such that vertices after the last midpoint are still assigned

    // assign 1D vertex value to all 3D vertices in vicinity
    for (size_t i = 0; i < outSize; i++) {
      auto vertexCoords = output()->vertices()[i].getCoords()[effectiveCoordinate];
      int  index        = 0;
      while (vertexCoords > axisMidpoints(index)) {
        index++;
        PRECICE_ASSERT(index < static_cast<int>(inSize));
      }
      outputValues((i * valueDimensions)) = inputValues(index);
    }
  } else {
    PRECICE_ASSERT(_type == COLLECT);
    /*
      3D vertices are projected onto the 1D axis and the data is then mapped
      to (and averaged at) the nearest 1D vertex in projection space.
    */
    PRECICE_ASSERT(outputValues.size() == static_cast<int>(output()->vertices().size()), outputValues.size(), valueDimensions, output()->vertices().size());
    size_t const inSize  = input()->vertices().size();
    size_t const outSize = output()->vertices().size();

    // determine principle axis midpoints as borders to assign 3D vertices to their respective 1D vertex
    Eigen::VectorXd axisMidpoints(outSize);
    for (size_t i = 0; i < (outSize - 1); i++) {
      auto axisPositionCurrent = output()->vertices()[i].getCoords()[effectiveCoordinate];
      auto axisPositionNext    = output()->vertices()[i + 1].getCoords()[effectiveCoordinate];
      axisMidpoints(i)         = (axisPositionCurrent + axisPositionNext) / 2;
    }
    axisMidpoints(outSize - 1) = std::numeric_limits<double>::max(); // large number, such that vertices after the last midpoint are still assigned

    Eigen::VectorXd counter(outSize); // counts number of vertices in between midpoints for averaging
    for (size_t i = 0; i < outSize; i++) {
      outputValues((i * valueDimensions)) = 0;
      counter(i)                          = 0;
    }

    // assign the 1D vertex the average of all 3D vertex values in vicinity
    for (size_t i = 0; i < inSize; i++) {
      auto vertexCoords = input()->vertices()[i].getCoords()[effectiveCoordinate];
      int  index        = 0;
      while (vertexCoords > axisMidpoints(index)) {
        index++;
        PRECICE_ASSERT(index < static_cast<int>(outSize));
      }
      outputValues(index) += inputValues(i * valueDimensions);
      counter(index) += 1;
    }
    for (size_t i = 0; i < outSize; i++) {
      outputValues(i) = outputValues(i) / counter(i);
    }
  }
}

void RadialGeoMultiscaleMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();

  computeMapping();

  if (getConstraint() == CONSISTENT) {
    PRECICE_ASSERT(_type == SPREAD, "Not yet implemented");

    // tag all vertices of the 1D mesh
    size_t const inSize = input()->vertices().size();
    for (size_t i = 0; i < inSize; i++) {
      input()->vertices()[i].tag();
    }

  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_ASSERT(false, "Not yet implemented");
  }

  clear();
}

void RadialGeoMultiscaleMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // no operation needed here for the moment
}

// TODO: needed for porting to develop
std::string RadialGeoMultiscaleMapping::getName() const
{
  return "radial-geomultiscale";
}

} // namespace precice::mapping
