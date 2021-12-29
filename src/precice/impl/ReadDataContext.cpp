#include "ReadDataContext.hpp"

namespace precice {
namespace impl {
ReadDataContext::ReadDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh)
    : DataContext(data, mesh)
{
}

void ReadDataContext::configureForReadMapping(MappingContext mappingContext, MeshContext fromMeshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(fromMeshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData fromData = fromMeshContext.mesh->data(getDataName());
  PRECICE_ASSERT(fromData != _providedData);
  this->setMapping(mappingContext, fromData, _providedData);
  PRECICE_ASSERT(hasReadMapping());
}

} // namespace impl
} // namespace precice