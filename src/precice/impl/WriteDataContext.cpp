#include "WriteDataContext.hpp"

namespace precice::impl {

logging::Logger WriteDataContext::_log{"impl::WriteDataContext"};

WriteDataContext::WriteDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh)
    : DataContext(data, mesh)
{
}

mesh::PtrData WriteDataContext::providedData()
{
  PRECICE_ASSERT(_providedData);
  return _providedData;
}

void WriteDataContext::writeValues(::precice::span<const VertexID> vertices, ::precice::span<const double> values)
{
  Eigen::Map<const Eigen::MatrixXd> inputData(values.data(), getDataDimensions(), vertices.size());
  Eigen::Map<Eigen::MatrixXd>       localData(_providedData->values().data(), getDataDimensions(), getMesh().vertices().size());

  for (int i = 0; i < vertices.size(); ++i) {
    const auto vid = vertices[i];
    PRECICE_CHECK(getMesh().isValidVertexID(vid),
                  "Cannot write data \"{}\" to invalid Vertex ID ({}) of mesh \"{}\". Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  getDataName(), vid, getMeshName());
    localData.col(vid) = inputData.col(i);
  }
}

void WriteDataContext::writeGradientValues(::precice::span<const VertexID> vertices, ::precice::span<const double> gradients)
{
  const auto                        gradientComponents = getMesh().getDimensions() * getDataDimensions();
  Eigen::Map<const Eigen::MatrixXd> inputGradients(gradients.data(), gradientComponents, vertices.size());
  Eigen::Map<Eigen::MatrixXd>       localGradients(_providedData->gradientValues().data(), gradientComponents, getMesh().vertices().size());

  for (int i = 0; i < vertices.size(); ++i) {
    const auto vid = vertices[i];
    PRECICE_CHECK(getMesh().isValidVertexID(vid),
                  "Cannot write gradient for data \"{}\" to invalid Vertex ID ({}) of mesh \"{}\". Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  getDataName(), vid, getMeshName());
    localGradients.col(vid) = inputGradients.col(i);
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
