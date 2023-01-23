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

  PRECICE_ASSERT(valueDimensions == output()->data(outputDataID)->getDimensions(),
                 valueDimensions, output()->data(outputDataID)->getDimensions());
  PRECICE_ASSERT(inputValues.size() / valueDimensions == static_cast<int>(input()->vertices().size()),
                 inputValues.size(), valueDimensions, input()->vertices().size());
  PRECICE_ASSERT(outputValues.size() / valueDimensions == static_cast<int>(output()->vertices().size()),
                 outputValues.size(), valueDimensions, output()->vertices().size());

  PRECICE_DEBUG("Map consistent");
  if (_type == SPREAD) {
    PRECICE_ASSERT(inputValues.size() == valueDimensions);
    if (_axis == X) {
      size_t const    inSize  = input()->vertices().size();
      size_t const    outSize = output()->vertices().size();
      Eigen::VectorXd axisMidpoints(inSize - 1);
      for (size_t i = 0; i < (inSize - 1); i++) {
        axisMidpoints(i) = (input()->vertices()[i].getCoords()[0] + input()->vertices()[i + 1].getCoords()[0]) / 2;
      }
      size_t i = 0;
      size_t j = 0;
      while (output()->vertices()[i].getCoords()[0] <= axisMidpoints(j)) {
        PRECICE_ASSERT(i <= outSize);
        PRECICE_ASSERT(j <= inSize);
        outputValues((i * valueDimensions) + 0) = inputValues((j * valueDimensions) + 0);
        i++;
        if (output()->vertices()[i + 1].getCoords()[0] > axisMidpoints(j)) {
          j++;
        }
      }
    } else if (_axis == Y) {
      size_t const    inSize  = input()->vertices().size();
      size_t const    outSize = output()->vertices().size();
      Eigen::VectorXd axisMidpoints(inSize - 1);
      for (size_t i = 0; i < (inSize - 1); i++) {
        axisMidpoints(i) = (input()->vertices()[i].getCoords()[1] + input()->vertices()[i + 1].getCoords()[1]) / 2;
      }
      size_t i = 0;
      size_t j = 0;
      while (output()->vertices()[i].getCoords()[1] <= axisMidpoints(j)) {
        PRECICE_ASSERT(i <= outSize);
        PRECICE_ASSERT(j <= inSize);
        outputValues((i * valueDimensions) + 1) = inputValues((j * valueDimensions) + 1);
        i++;
        if (output()->vertices()[i + 1].getCoords()[1] > axisMidpoints(j)) {
          j++;
        }
      }
    } else {
      PRECICE_ASSERT(_axis == Z);
      size_t const    inSize  = input()->vertices().size();
      size_t const    outSize = output()->vertices().size();
      Eigen::VectorXd axisMidpoints(inSize - 1);
      for (size_t i = 0; i < (inSize - 1); i++) {
        axisMidpoints(i) = (input()->vertices()[i].getCoords()[2] + input()->vertices()[i + 1].getCoords()[2]) / 2;
      }
      size_t i = 0;
      size_t j = 0;
      while (output()->vertices()[i].getCoords()[2] <= axisMidpoints(j)) {
        PRECICE_ASSERT(i <= outSize);
        PRECICE_ASSERT(j <= inSize);
        outputValues((i * valueDimensions) + 2) = inputValues((j * valueDimensions) + 2);
        i++;
        if (output()->vertices()[i + 1].getCoords()[2] > axisMidpoints(j)) {
          j++;
        }
      }
    }
  } else {
    PRECICE_ASSERT(_type == COLLECT);
    PRECICE_ASSERT(outputValues.size() == valueDimensions);
    //PRECICE_ASSERT(output()->vertices().size() == 1);
    //for (int dim = 0; dim < valueDimensions; dim++) {
    //  outputValues(dim) = 0.0;
    //}
    //size_t const inSize = input()->vertices().size();
    //for (size_t i = 0; i < inSize; i++) {
    //  for (int dim = 0; dim < valueDimensions; dim++) {
    //    outputValues(dim) += inputValues((i * valueDimensions) + dim) / inSize;
    //  }
    //}
    if (_axis == X) {
      size_t const    inSize  = input()->vertices().size();
      size_t const    outSize = output()->vertices().size();
      Eigen::VectorXd axisMidpoints(inSize - 1);
      for (size_t i = 0; i < (outSize - 1); i++) {
        axisMidpoints(i) = (output()->vertices()[i].getCoords()[0] + output()->vertices()[i + 1].getCoords()[0]) / 2;
      }
      size_t i       = 0;
      size_t j       = 0;
      size_t counter = 0;
      for (size_t index = 0; index < outSize; index++) {
        outputValues((index * valueDimensions) + 0) = 0;
      }
      while (input()->vertices()[i].getCoords()[0] <= axisMidpoints(j)) {
        PRECICE_ASSERT(i <= inSize);
        PRECICE_ASSERT(j <= outSize);
        outputValues((j * valueDimensions) + 0) += inputValues((i * valueDimensions) + 0);
        i++;
        counter++;
        if (input()->vertices()[i + 1].getCoords()[0] > axisMidpoints(j)) {
          j++;
          outputValues((j * valueDimensions) + 0) = outputValues((j * valueDimensions) + 0) / counter;
          counter                                 = 0;
        }
      }
    } else if (_axis == Y) {
      size_t const    inSize  = input()->vertices().size();
      size_t const    outSize = output()->vertices().size();
      Eigen::VectorXd axisMidpoints(inSize - 1);
      for (size_t i = 0; i < (outSize - 1); i++) {
        axisMidpoints(i) = (output()->vertices()[i].getCoords()[1] + output()->vertices()[i + 1].getCoords()[1]) / 2;
      }
      size_t i       = 0;
      size_t j       = 0;
      size_t counter = 0;
      for (size_t index = 0; index < outSize; index++) {
        outputValues((index * valueDimensions) + 1) = 0;
      }
      while (input()->vertices()[i].getCoords()[1] <= axisMidpoints(j)) {
        PRECICE_ASSERT(i <= inSize);
        PRECICE_ASSERT(j <= outSize);
        outputValues((j * valueDimensions) + 1) += inputValues((i * valueDimensions) + 1);
        i++;
        counter++;
        if (input()->vertices()[i + 1].getCoords()[1] > axisMidpoints(j)) {
          j++;
          outputValues((j * valueDimensions) + 1) = outputValues((j * valueDimensions) + 1) / counter;
          counter                                 = 0;
        }
      }
    } else {
      PRECICE_ASSERT(_axis == Z);
      size_t const    inSize  = input()->vertices().size();
      size_t const    outSize = output()->vertices().size();
      Eigen::VectorXd axisMidpoints(inSize - 1);
      for (size_t i = 0; i < (outSize - 1); i++) {
        axisMidpoints(i) = (output()->vertices()[i].getCoords()[1] + output()->vertices()[i + 1].getCoords()[1]) / 2;
      }
      size_t i       = 0;
      size_t j       = 0;
      size_t counter = 0;
      for (size_t index = 0; index < outSize; index++) {
        outputValues((index * valueDimensions) + 1) = 0;
      }
      while (input()->vertices()[i].getCoords()[1] <= axisMidpoints(j)) {
        PRECICE_ASSERT(i <= inSize);
        PRECICE_ASSERT(j <= outSize);
        outputValues((j * valueDimensions) + 1) += inputValues((i * valueDimensions) + 1);
        i++;
        counter++;
        if (input()->vertices()[i + 1].getCoords()[1] > axisMidpoints(j)) {
          j++;
          outputValues((j * valueDimensions) + 1) = outputValues((j * valueDimensions) + 1) / counter;
          counter                                 = 0;
        }
      }
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

} // namespace mapping
} // namespace precice