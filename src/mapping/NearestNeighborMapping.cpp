#include "NearestNeighborMapping.hpp"

#include <Eigen/Core>
#include <boost/container/flat_set.hpp>
#include <functional>
#include "logging/LogMacros.hpp"
#include "profiling/Event.hpp"
#include "profiling/EventUtils.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::mapping {

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

void NearestNeighborMapping::evaluateMappingDataCacheAt(::precice::span<const double> coordinates, const MappingDataCache &cache, ::precice::span<double> values)
{
  auto  searchSpace = input();
  auto &index       = searchSpace->index();
  auto  dim         = getDimensions();

  // Set up of output arrays
  Eigen::Map<const Eigen::MatrixXd> localData(cache.inData.data(), cache.getDataDimensions(), cache.inData.size() / cache.getDataDimensions());
  Eigen::Map<Eigen::MatrixXd>       outputData(values.data(), cache.getDataDimensions(), values.size());

  const size_t verticesSize = coordinates.size() / dim;
  for (size_t i = 0; i < verticesSize; ++i) {
    Eigen::Map<const Eigen::VectorXd> localCoords(coordinates.data() + i * dim, dim);
    const auto &                      matchedVertex = index.getClosestVertex(localCoords);
    outputData.col(i)                               = localData.col(matchedVertex.index);
  }
}

void NearestNeighborMapping::writeConservativeAt(::precice::span<const double> coordinates, Eigen::Map<const Eigen::MatrixXd> &source, Eigen::Map<Eigen::MatrixXd> &target)
{
  auto  searchSpace = output();
  auto &index       = searchSpace->index();
  auto  dim         = getDimensions();

  const size_t verticesSize = coordinates.size() / dim;
  for (size_t i = 0; i < verticesSize; ++i) {
    Eigen::Map<const Eigen::VectorXd> localCoords(coordinates.data() + i * dim, dim);
    const auto &                      matchedVertex = index.getClosestVertex(localCoords);
    target.col(matchedVertex.index) += source.col(i);
  }
}

void NearestNeighborMapping::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map." + mappingNameShort + ".mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_DEBUG("Map conservative using {}", getName());

  const Eigen::VectorXd &inputValues  = inData.values;
  Eigen::VectorXd &      outputValues = outData;

  // Data dimensions (for scalar = 1, for vectors > 1)
  const size_t inSize          = input()->nVertices();
  const int    valueDimensions = inData.dataDims;

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

void NearestNeighborMapping::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map." + mappingNameShort + ".mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_DEBUG("Map {} using {}", (hasConstraint(CONSISTENT) ? "consistent" : "scaled-consistent"), getName());

  const Eigen::VectorXd &inputValues  = inData.values;
  Eigen::VectorXd &      outputValues = outData;

  // Data dimensions (for scalar = 1, for vectors > 1)
  const size_t outSize         = output()->nVertices();
  const int    valueDimensions = inData.dataDims;

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

} // namespace precice::mapping
