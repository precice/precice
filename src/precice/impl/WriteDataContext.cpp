#include "WriteDataContext.hpp"

namespace precice {
namespace impl {

logging::Logger WriteDataContext::_log{"impl::WriteDataContext"};

WriteDataContext::WriteDataContext(mesh::PtrData data, mesh::PtrMesh mesh) : DataContext(data, mesh) {}

mesh::PtrData WriteDataContext::providedData()
{
  PRECICE_ASSERT(_providedData);
  return _providedData;
}

void WriteDataContext::appendMappingConfiguration(const MappingContext &mappingContext, const MeshContext &meshContext)
{
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData,
                 "Data the write mapping is mapping to needs to be different from _providedData");
  appendMapping(mappingContext, _providedData, data);
  PRECICE_ASSERT(hasWriteMapping());
}

} // namespace impl
} // namespace precice
