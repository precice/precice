
#include "NearestNeighborGradientMapping.hpp"

#include <Eigen/Core>
#include <boost/container/flat_set.hpp>
#include <functional>
#include <iostream>
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Event.hpp"
#include "utils/assertion.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

NearestNeighborGradientMapping::NearestNeighborGradientMapping(
    Constraint constraint,
    int        dimensions)
    : NearestNeighborBaseMapping(constraint, dimensions, true, "NearestNeighborGradientMapping", "nng")
{
  if (hasConstraint(SCALEDCONSISTENT)) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    setInputRequirement(Mapping::MeshRequirement::GRADIENT);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
}

void NearestNeighborGradientMapping::map(
    int inputDataID,
    int outputDataID)
{
  // TODO: CHANGE TO DEBUG (BETTER DEBUG INFOS)
  PRECICE_TRACE(inputDataID, outputDataID);

  precice::utils::Event e("map." + MAPPING_NAME_SHORT + ".mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  int valueDimensions = input()->data(inputDataID)->getDimensions(); // Data dimensions (bei scalar = 1, bei vectors > 1)

  const Eigen::VectorXd &inputValues  = input()->data(inputDataID)->values();
  Eigen::VectorXd &      outputValues = output()->data(outputDataID)->values();

  /// Check if input has gradient data, else send Error
  /// TODO: NN WAS REQUESTED AND HOW TO FIX THIS
  //TODO: INPUT SHOULD BE EMPTY OR NOT ? DIVIDE THE TWO CHECKS !!
  // TODO: MAYBE USE A PRECICE_CHECK ?
  if (!input()->vertices().empty() && !input()->data(inputDataID)->hasGradient()) {
    PRECICE_ERROR("Mesh \"{}\" does not contain gradient data. ",
                  input()->getName());
  }

  const Eigen::MatrixXd &gradientValues = input()->data(inputDataID)->gradientValues();

  //assign(outputValues) = 0.0;

  //TODO: REMOVE THIS ASSERT ?
  PRECICE_ASSERT(valueDimensions == output()->data(outputDataID)->getDimensions(),
                 valueDimensions, output()->data(outputDataID)->getDimensions())
  PRECICE_ASSERT(inputValues.size() / valueDimensions == (int) input()->vertices().size(),
                 inputValues.size(), valueDimensions, input()->vertices().size());
  PRECICE_ASSERT(outputValues.size() / valueDimensions == (int) output()->vertices().size(),
                 outputValues.size(), valueDimensions, output()->vertices().size());

  if (hasConstraint(CONSERVATIVE)) {
    PRECICE_DEBUG("Map conservative");
    const size_t inSize = input()->vertices().size();

    for (size_t i = 0; i < inSize; i++) {
      const int outputIndex = _vertexIndices[i] * valueDimensions;

      for (int dim = 0; dim < valueDimensions; dim++) {

        int mapOutputIndex = outputIndex + dim;
        int mapInputIndex  = (i * valueDimensions) + dim;

        outputValues(mapOutputIndex) += inputValues(mapInputIndex) + _distancesMatched[i].transpose() * gradientValues.col(mapInputIndex);
      }
    }
  } else {
    PRECICE_DEBUG((hasConstraint(CONSISTENT) ? "Map consistent" : "Map scaled-consistent"));
    size_t const outSize = output()->vertices().size();

    for (size_t i = 0; i < outSize; i++) {
      int inputIndex = _vertexIndices[i] * valueDimensions;

      for (int dim = 0; dim < valueDimensions; dim++) {

        int mapOutputIndex = (i * valueDimensions) + dim;
        int mapInputIndex  = inputIndex + dim;

        outputValues(mapOutputIndex) = inputValues(mapInputIndex) + _distancesMatched[i].transpose() * gradientValues.col(mapInputIndex);
      }
    }
    if (hasConstraint(SCALEDCONSISTENT)) {
      scaleConsistentMapping(inputDataID, outputDataID);
    }

    PRECICE_DEBUG("Mapped values (with gradient) = {}", utils::previewRange(3, outputValues));
  }
}

} // namespace mapping
} // namespace precice
