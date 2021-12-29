#include "WriteDataContext.hpp"

namespace precice {
namespace impl {
WriteDataContext::WriteDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh)
    : DataContext(data, mesh)
{
}

void WriteDataContext::configureForWriteMapping(MappingContext mappingContext, MeshContext toMeshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(toMeshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData toData = toMeshContext.mesh->data(getDataName());
  PRECICE_ASSERT(toData != _providedData);
  this->setMapping(mappingContext, _providedData, toData);
  PRECICE_ASSERT(hasWriteMapping());
}

} // namespace impl
} // namespace precice