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

void NearestNeighborMapping::mapConsistentAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> values)
{
  precice::profiling::Event e("map.nn.mapConsistentAt.From" + input()->getName());
  auto &                    index = input()->index();

  // Set up of output arrays
  Eigen::Map<const Eigen::MatrixXd> localData(cache.inData.data(), cache.getDataDimensions(), cache.inData.size() / cache.getDataDimensions());
  for (Eigen::Index i = 0; i < coordinates.cols(); ++i) {
    values.col(i) = localData.col(index.getClosestVertex(coordinates.col(i)).index);
  }
}

void NearestNeighborMapping::mapConservativeAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, impl::MappingDataCache &cache, const Eigen::Ref<const Eigen::MatrixXd> &source, Eigen::Ref<Eigen::MatrixXd> target)
{
  precice::profiling::Event e("map.nn.mapConservativeAt.From" + input()->getName());
  auto &                    index = output()->index();
  for (Eigen::Index i = 0; i < coordinates.cols(); ++i) {
    target.col(index.getClosestVertex(coordinates.col(i)).index) += source.col(i);
  }
}

void NearestNeighborMapping::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map." + mappingNameShort + ".mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_DEBUG("Map conservative using {}", getName());

  // scalar = 1, vector > 1
  const int valueDimensions = inData.dataDims;
  // Create Eigen::Map to interpret data as a matrix (enables col-wise operations)
  Eigen::Map<const Eigen::MatrixXd> inMap(inData.values.data(), valueDimensions, input()->nVertices());
  Eigen::Map<Eigen::MatrixXd>       outMap(outData.data(), valueDimensions, output()->nVertices());

  // Apply mapping
  for (Eigen::Index i = 0; i < inMap.cols(); ++i) {
    outMap.col(_vertexIndices[i]) += inMap.col(i);
  }

  PRECICE_DEBUG("Mapped values = {}", utils::previewRange(3, outData));
}

void NearestNeighborMapping::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map." + mappingNameShort + ".mapData.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);
  PRECICE_DEBUG("Map {} using {}", (hasConstraint(CONSISTENT) ? "consistent" : "scaled-consistent"), getName());

  const int valueDimensions = inData.dataDims;
  // Map the raw pointers as (valueDimensions x vertices) in column-major Eigen style
  Eigen::Map<const Eigen::MatrixXd> inputMap(inData.values.data(), valueDimensions, input()->nVertices());
  Eigen::Map<Eigen::MatrixXd>       outputMap(outData.data(), valueDimensions, output()->nVertices());

  for (Eigen::Index i = 0; i < outputMap.cols(); ++i) {
    // _vertexIndices[i] is the solution of our mapping
    outputMap.col(i) = inputMap.col(_vertexIndices[i]);
  }
  PRECICE_DEBUG("Mapped values = {}", utils::previewRange(3, outData));
}

std::string NearestNeighborMapping::getName() const
{
  return "nearest-neighbor";
}

} // namespace precice::mapping
