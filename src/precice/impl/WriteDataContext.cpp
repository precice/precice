#include "WriteDataContext.hpp"

namespace precice {
namespace impl {
WriteDataContext::WriteDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh)
    : DataContext(data, mesh)
{
}

void WriteDataContext::configureMapping(MappingContext mappingContext, MeshContext meshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the write mapping is mapping to needs to be different from _providedData");
  setMapping(mappingContext, _providedData, data);
  PRECICE_ASSERT(hasWriteMapping());
}

} // namespace impl
} // namespace precice