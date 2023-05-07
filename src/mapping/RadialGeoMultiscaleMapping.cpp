#include "RadialGeoMultiscaleMapping.hpp"
#include "mesh/Mesh.hpp"

namespace precice {

namespace mapping {

RadialGeoMultiscaleMapping::RadialGeoMultiscaleMapping(
    Constraint     constraint,
    int            dimensions,
    MultiscaleType type,
    RadialAxis     axis)
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
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
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
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
  PRECICE_ASSERT(false, "Not yet implemented");
  PRECICE_DEBUG("Map conservative");
}

void RadialGeoMultiscaleMapping::mapConsistent(DataID inputDataID, DataID outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);

  const Eigen::VectorXd &inputValues     = input()->data(inputDataID)->values();
  int                    valueDimensions = input()->data(inputDataID)->getDimensions();
  Eigen::VectorXd &      outputValues    = output()->data(outputDataID)->values();

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

  PRECICE_ASSERT((inputValues.size() / valueDimensions == static_cast<int>(input()->vertices().size())),
                 inputValues.size(), valueDimensions, input()->vertices().size());
  PRECICE_ASSERT((outputValues.size() / valueDimensions == static_cast<int>(output()->vertices().size())),
                 outputValues.size(), valueDimensions, output()->vertices().size());
  PRECICE_ASSERT(valueDimensions == 1);

  PRECICE_DEBUG("Map consistent");
  if (_type == SPREAD) {
    size_t const inSize  = input()->vertices().size();
    size_t const outSize = output()->vertices().size();

    Eigen::VectorXd axisMidpoints(inSize);
    for (size_t i = 0; i < (inSize - 1); i++) {
      auto axisPositionCurrent = input()->vertices()[i].getCoords()[coord];
      auto axisPositionNext    = input()->vertices()[i + 1].getCoords()[coord];
      axisMidpoints(i)         = (axisPositionCurrent + axisPositionNext) / 2;
    }
    axisMidpoints(inSize - 1) = 1e100; // arbitrary large number, such that vertices after the last midpoint are still assigned
    for (size_t i = 0; i < outSize; i++) {
      auto vertexCoords = output()->vertices()[i].getCoords()[coord];
      int  index        = 0;
      while (vertexCoords > axisMidpoints(index)) {
        index++;
        PRECICE_ASSERT(index < static_cast<int>(inSize));
      }
      outputValues((i * valueDimensions)) = inputValues(index);
    }
  } else {
    PRECICE_ASSERT(_type == COLLECT);
    PRECICE_ASSERT(outputValues.size() == static_cast<int>(output()->vertices().size()), outputValues.size(), valueDimensions, output()->vertices().size());
    size_t const inSize  = input()->vertices().size();
    size_t const outSize = output()->vertices().size();

    Eigen::VectorXd axisMidpoints(outSize);
    for (size_t i = 0; i < (outSize - 1); i++) {
      auto axisPositionCurrent = output()->vertices()[i].getCoords()[coord];
      auto axisPositionNext    = output()->vertices()[i + 1].getCoords()[coord];
      axisMidpoints(i)         = (axisPositionCurrent + axisPositionNext) / 2;
    }
    axisMidpoints(outSize - 1) = 1e100; // arbitrary large number, such that vertices after the last midpoint are still assigned

    Eigen::VectorXd counter(outSize); // counts number of vertices inbetween midpoints for averaging
    for (size_t i = 0; i < outSize; i++) {
      outputValues((i * valueDimensions)) = 0;
      counter(i)                          = 0;
    }

    for (size_t i = 0; i < inSize; i++) {
      auto vertexCoords = input()->vertices()[i].getCoords()[coord];
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
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE, getConstraint());
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
// std::string AxialGeoMultiscaleMapping::getName() const
// {
//   return "radial-geomultiscale";
// }

} // namespace mapping
} // namespace precice