#include "ReadDataContext.hpp"

namespace precice::impl {

logging::Logger ReadDataContext::_log{"impl::ReadDataContext"};

ReadDataContext::ReadDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh)
    : DataContext(data, mesh)
{
}

void ReadDataContext::appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext)
{
  PRECICE_ASSERT(!hasReadMapping(), "The read data context must be unique. Otherwise we would have an ambiguous read data operation on the user side.");
  // The read mapping must be unique, but having read and write in the same context is not possible either
  PRECICE_ASSERT(_mappingContexts.empty());
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the read mapping is mapping from needs to be different from _providedData");
  mappingContext.fromData = data;
  mappingContext.toData   = _providedData;
  appendMapping(mappingContext);
  PRECICE_ASSERT(hasReadMapping());
}

bool ReadDataContext::hasSamples() const
{
  return _providedData->hasSamples();
}

void ReadDataContext::readValues(::precice::span<const VertexID> vertices, double readTime, ::precice::span<double> values) const
{
  Eigen::Map<Eigen::MatrixXd> outputData(values.data(), getDataDimensions(), values.size());
  auto                        sampleResult = _providedData->sampleAtTime(readTime);
  auto                        localData    = sampleResult.values().reshaped(getDataDimensions(), getMeshVertexCount());
  for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
    outputData.col(i) = localData.col(vertices[i]);
  }
}

void ReadDataContext::mapAndReadValues(::precice::span<const double> coordinates, double readTime, ::precice::span<double> values)
{
  PRECICE_TRACE(getMeshName(), getDataName(), coordinates.size(), values.size(), readTime);
  PRECICE_ASSERT(mappingCache);
  PRECICE_ASSERT(justInTimeMapping);

  if (!mappingCache->hasDataAtTimeStamp(readTime)) {
    // Sample waveform relaxation
    Eigen::VectorXd sample{_providedData->sampleAtTime(readTime)};
    justInTimeMapping->updateMappingDataCache(*mappingCache.get(), sample);
    mappingCache->setTimeStamp(readTime);
  }
  Eigen::Map<const Eigen::MatrixXd> coords(coordinates.data(), getSpatialDimensions(), coordinates.size() / getSpatialDimensions());
  Eigen::Map<Eigen::MatrixXd>       target(values.data(), getDataDimensions(), values.size() / getDataDimensions());

  // Function, which fills the values using the coordinates and the cache
  justInTimeMapping->mapConsistentAt(coords, *mappingCache.get(), target);
}

int ReadDataContext::getWaveformDegree() const
{
  return _providedData->getWaveformDegree();
}

void ReadDataContext::clearToDataFor(const cplscheme::ImplicitData &from)
{
  PRECICE_TRACE(getMeshName(), getDataName());
  PRECICE_ASSERT(hasMapping());
  for (auto &context : _mappingContexts) {
    auto id = context.fromData->getID();
    if (from.contains(id)) {
      if (from.toKeep(id)) {
        context.toData->timeStepsStorage().clearExceptLast();
      } else {
        context.toData->timeStepsStorage().clear();
      }
    }
  }
}

void ReadDataContext::trimToDataAfterFor(const cplscheme::ImplicitData &from, double t)
{
  PRECICE_TRACE(getMeshName(), getDataName(), t);
  PRECICE_ASSERT(hasMapping());
  for (auto &context : _mappingContexts) {
    if (from.contains(context.fromData->getID())) {
      context.toData->timeStepsStorage().trimAfter(t);
    }
  }
}

} // namespace precice::impl
