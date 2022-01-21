#include "ReadDataContext.hpp"

namespace precice {
namespace impl {
ReadDataContext::ReadDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh)
    : DataContext(data, mesh)
{
}

void ReadDataContext::configureMapping(MappingContext mappingContext, MeshContext meshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the read mapping is mapping from needs to be different from _providedData");
  setMapping(mappingContext, data, _providedData);
  PRECICE_ASSERT(hasReadMapping());
}

} // namespace impl
} // namespace precice