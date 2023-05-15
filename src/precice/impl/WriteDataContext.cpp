#include "WriteDataContext.hpp"

namespace precice::impl {

logging::Logger WriteDataContext::_log{"impl::WriteDataContext"};

WriteDataContext::WriteDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh)
    : DataContext(data, mesh)
{
}

void WriteDataContext::writeValues(::precice::span<const VertexID> vertices, ::precice::span<const double> values)
{
  Eigen::Map<const Eigen::MatrixXd> inputData(values.data(), getDataDimensions(), vertices.size());
  Eigen::Map<Eigen::MatrixXd>       localData(_providedData->values().data(), getDataDimensions(), getMeshVertexCount());

  for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
    localData.col(vertices[i]) = inputData.col(i);
  }
}

void WriteDataContext::writeGradients(::precice::span<const VertexID> vertices, ::precice::span<const double> gradients)
{
  const auto                        gradientComponents = getSpatialDimensions() * getDataDimensions();
  Eigen::Map<const Eigen::MatrixXd> inputGradients(gradients.data(), gradientComponents, vertices.size());
  Eigen::Map<Eigen::MatrixXd>       localGradients(_providedData->gradients().data(), gradientComponents, getMeshVertexCount());

  for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
    localGradients.col(vertices[i]) = inputGradients.col(i);
  }
}

void WriteDataContext::appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext)
{
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the write mapping is mapping to needs to be different from _providedData");
  mappingContext.fromData = _providedData;
  mappingContext.toData   = data;
  appendMapping(mappingContext);
  PRECICE_ASSERT(hasWriteMapping());
}

} // namespace precice::impl
