
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
  PRECICE_CHECK(!hasConstraint(CONSERVATIVE), "Nearest-neighbor-gradient mapping is not implemented using a \"conservative\" constraint. Please select constraint=\""consistent\" or a different mapping method.");

  if (hasConstraint(SCALEDCONSISTENT)) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
}

void NearestNeighborGradientMapping::onMappingComputed(mesh::PtrMesh origins, mesh::PtrMesh searchSpace)
{

  // Initialize the offsets list
  _offsetsMatched.resize(_vertexIndices.size());

  // Calculate offsets
  for (size_t i = 0; i < _vertexIndices.size(); ++i) {

    const auto &matchedVertexCoords = searchSpace.get()->vertices()[_vertexIndices[i]].getCoords();
    const auto &sourceVertexCoords  = origins->vertices()[i].getCoords();

    // We calculate the distances uniformly for consistent mapping constraint as the difference (output - input)
    _offsetsMatched[i] = sourceVertexCoords - matchedVertexCoords;
  }
};

void NearestNeighborGradientMapping::map(
    int inputDataID,
    int outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);

  precice::utils::Event e("map." + mappingNameShort + ".mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  const int valueDimensions = input()->data(inputDataID)->getDimensions(); // Data dimensions (for scalar = 1, for vectors > 1)

  const Eigen::VectorXd &inputValues  = input()->data(inputDataID)->values();
  Eigen::VectorXd &      outputValues = output()->data(outputDataID)->values();

  /// Check if input has gradient data, else send Error
  if (input()->vertices().empty()) {
    PRECICE_WARN("The mesh doesn't contain any vertices.");
  }

  PRECICE_CHECK(input()->data(inputDataID)->hasGradient(), "Mesh \"{}\" does not contain gradient data. Using Nearest Neighbor Gradient requires gradient data for each vertices.",
                "Check if hasGradient flag in the Data object was successfully initialized.",
                input()->getName());

  const Eigen::MatrixXd &gradientValues = input()->data(inputDataID)->gradientValues();

  PRECICE_ASSERT(inputValues.size() / valueDimensions == static_cast<int>(input()->vertices().size()),
                 inputValues.size(), valueDimensions, input()->vertices().size());
  PRECICE_ASSERT(outputValues.size() / valueDimensions == static_cast<int>(output()->vertices().size()),
                 outputValues.size(), valueDimensions, output()->vertices().size());

  //Consistent mapping
  PRECICE_DEBUG((hasConstraint(CONSISTENT) ? "Map consistent" : "Map scaled-consistent"));
  size_t const outSize = output()->vertices().size();

  for (size_t i = 0; i < outSize; i++) {
    int inputIndex = _vertexIndices[i] * valueDimensions;

    for (int dim = 0; dim < valueDimensions; dim++) {

      const int mapOutputIndex = (i * valueDimensions) + dim;
      const int mapInputIndex  = inputIndex + dim;

      outputValues(mapOutputIndex) = inputValues(mapInputIndex) + _offsetsMatched[i].transpose() * gradientValues.col(mapInputIndex);
    }
  }

  if (hasConstraint(SCALEDCONSISTENT)) {
    scaleConsistentMapping(inputDataID, outputDataID);
  }

  PRECICE_DEBUG("Mapped values (with gradient) = {}", utils::previewRange(3, outputValues));
}

} // namespace mapping
} // namespace precice