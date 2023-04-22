#include "NearestNeighborMapping.hpp"

#include <Eigen/Core>
#include <boost/container/flat_set.hpp>
#include <functional>
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Event.hpp"
#include "utils/EventUtils.hpp"
#include "utils/assertion.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

NearestNeighborMapping::NearestNeighborMapping(
    Constraint constraint,
    int        dimensions)
    : NearestNeighborBaseMapping(constraint, dimensions, false, "NearestNeighborMapping", "nn")
{
  if (isScaledConsistent()) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
}

void NearestNeighborMapping::mapConservative(DataID inputDataID, DataID outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);
  precice::utils::Event e("map." + mappingNameShort + ".mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
  PRECICE_DEBUG("Map conservative");

  const Eigen::VectorXd &inputValues  = input()->data(inputDataID)->values();
  Eigen::VectorXd &      outputValues = output()->data(outputDataID)->values();

  // Data dimensions (for scalar = 1, for vectors > 1)
  const int    valueDimensions = input()->data(inputDataID)->getDimensions();
  const size_t inSize          = input()->vertices().size();

  for (size_t i = 0; i < inSize; i++) {
    int const outputIndex = _vertexIndices[i] * valueDimensions;

    for (int dim = 0; dim < valueDimensions; dim++) {

      const int mapOutputIndex = outputIndex + dim;
      const int mapInputIndex  = (i * valueDimensions) + dim;

      outputValues(mapOutputIndex) += inputValues(mapInputIndex);
    }
  }
  PRECICE_DEBUG("Mapped values = {}", utils::previewRange(3, outputValues));
}

void NearestNeighborMapping::mapConsistent(DataID inputDataID, DataID outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);
  precice::utils::Event e("map." + mappingNameShort + ".mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
  PRECICE_DEBUG((hasConstraint(CONSISTENT) ? "Map consistent" : "Map scaled-consistent"));

  const Eigen::VectorXd &inputValues  = input()->data(inputDataID)->values();
  Eigen::VectorXd &      outputValues = output()->data(outputDataID)->values();

  // Data dimensions (for scalar = 1, for vectors > 1)
  const int    valueDimensions = input()->data(inputDataID)->getDimensions();
  const size_t outSize         = output()->vertices().size();

  for (size_t i = 0; i < outSize; i++) {
    int inputIndex = _vertexIndices[i] * valueDimensions;

    for (int dim = 0; dim < valueDimensions; dim++) {

      const int mapOutputIndex = (i * valueDimensions) + dim;
      const int mapInputIndex  = inputIndex + dim;

      outputValues(mapOutputIndex) = inputValues(mapInputIndex);
    }
  }
  PRECICE_DEBUG("Mapped values = {}", utils::previewRange(3, outputValues));
}

std::string NearestNeighborMapping::getName() const
{
  return "nearest-neighbor";
}
} // namespace mapping
} // namespace precice
