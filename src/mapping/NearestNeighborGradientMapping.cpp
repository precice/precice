#include "NearestNeighborGradientMapping.hpp"

#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <boost/container/flat_set.hpp>
#include <functional>
#include <iostream>
#include <strings.h>
#include "logging/LogMacros.hpp"
#include "profiling/Event.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::mapping {

NearestNeighborGradientMapping::NearestNeighborGradientMapping(
    Constraint constraint,
    int        dimensions)
    : NearestNeighborBaseMapping(constraint, dimensions, true, "NearestNeighborGradientMapping", "nng")
{
  PRECICE_ASSERT(!hasConstraint(CONSERVATIVE));

  if (isScaledConsistent()) {
    PRECICE_WARN("The scaled-consistent mapping hasn't been specifically tested with nearest-neighbor-gradient. Please avoid using it or choose another mapping method. ");
  }

  if (isScaledConsistent()) {
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
    // For consistent mapping: the source is the output vertex and the matched vertex is the input since we iterate over all outputs
    // and assign each exactly one vertex form the search space, which are our origins vertices.
    _offsetsMatched[i] = sourceVertexCoords - matchedVertexCoords;
  }
};

void NearestNeighborGradientMapping::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map." + mappingNameShort + ".mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);

  PRECICE_ASSERT(inData.values.size() == 0 || inData.gradients.size() != 0,
                 "Mesh \"{}\" does not contain gradient data. Using Nearest Neighbor Gradient mapping requires gradient data.",
                 input()->getName());

  /// Check if input has gradient data, else send Error
  if (input()->vertices().empty()) {
    PRECICE_WARN("The mesh doesn't contain any vertices.");
  }

  const int              valueDimensions = inData.dataDims;
  const Eigen::VectorXd &inputValues     = inData.values;
  Eigen::VectorXd &      outputValues    = outData;
  const Eigen::MatrixXd &gradients       = inData.gradients;

  // Consistent mapping
  PRECICE_DEBUG("Map {} using {}", (hasConstraint(CONSISTENT) ? "consistent" : "scaled-consistent"), getName());
  const size_t outSize = output()->vertices().size();

  for (size_t i = 0; i < outSize; i++) {
    int inputIndex = _vertexIndices[i] * valueDimensions;

    for (int dim = 0; dim < valueDimensions; dim++) {

      const int mapOutputIndex = (i * valueDimensions) + dim;
      const int mapInputIndex  = inputIndex + dim;

      outputValues(mapOutputIndex) = inputValues(mapInputIndex) + _offsetsMatched[i].transpose() * gradients.col(mapInputIndex);
    }
  }

  PRECICE_DEBUG("Mapped values (with gradient) = {}", utils::previewRange(3, outputValues));
}

void NearestNeighborGradientMapping::mapConservative(const time::Sample & /* inData */, Eigen::VectorXd & /* outData */)
{
  PRECICE_ASSERT(false, "Not implemented.");
}

std::string NearestNeighborGradientMapping::getName() const
{
  return "nearest-neighbor-gradient";
}

} // namespace precice::mapping
