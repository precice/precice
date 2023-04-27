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

void WriteDataContext::writeValues(std::vector<int> indices, const Eigen::Map<const Eigen::VectorXd> values)
{
  const auto vertexCount    = getDataSize() / getDataDimensions();
  auto &     valuesInternal = providedData()->values();
  for (int i = 0; i < indices.size(); i++) {
    const auto valueIndex = indices[i];
    PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                  "Cannot write data \"{}\" to invalid Vertex ID ({}). Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  getDataName(), valueIndex);
    const int offsetInternal = valueIndex * getDataDimensions();
    const int offset         = i * getDataDimensions();
    for (int dim = 0; dim < getDataDimensions(); dim++) {
      PRECICE_ASSERT(offset + dim < getDataSize(),
                     offset + dim, getDataSize());
      valuesInternal[offsetInternal + dim] = values[offset + dim];
    }
  }
}

void WriteDataContext::writeGradientValues(std::vector<int> indices, const Eigen::Map<const Eigen::MatrixXd> gradients)
{
  const auto vertexCount            = getDataSize() / getDataDimensions();
  auto &     gradientValuesInternal = providedData()->gradientValues();
  const int  stride                 = getDataDimensions();
  PRECICE_ASSERT(providedData() != nullptr);

  for (auto i = 0; i < indices.size(); i++) {
    const auto valueIndex = indices[i];
    PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                  "Cannot write gradient data \"{}\" to invalid Vertex ID ({}). Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  getDataName(), valueIndex);
    gradientValuesInternal.block(0, stride * valueIndex, getSpatialDimensions(), getDataDimensions()) = gradients.block(0, stride * i, getSpatialDimensions(), getDataDimensions());
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
