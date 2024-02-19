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
  PRECICE_ASSERT(!hasReadMapping(), "The read data context must be unique. Otherwise we would have an ambiguous read data operation on the user side.")
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the read mapping is mapping from needs to be different from _providedData");
  mappingContext.fromData = data;
  mappingContext.toData   = _providedData;
  appendMapping(mappingContext);
  PRECICE_ASSERT(hasReadMapping());
}

void ReadDataContext::readValues(::precice::span<const VertexID> vertices, double readTime, ::precice::span<double> values) const
{
  Eigen::Map<Eigen::MatrixXd>       outputData(values.data(), getDataDimensions(), values.size());
  const Eigen::MatrixXd             sample{_providedData->sampleAtTime(readTime)};
  Eigen::Map<const Eigen::MatrixXd> localData(sample.data(), getDataDimensions(), getMeshVertexCount());
  for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
    outputData.col(i) = localData.col(vertices[i]);
  }
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
