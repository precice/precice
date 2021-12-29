#include "WriteDataContext.hpp"

namespace precice {
namespace impl {
WriteDataContext::WriteDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh)
    : DataContext(data, mesh)
{
}

mesh::PtrData WriteDataContext::providedData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_providedData);
  return _providedData;
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

void WriteDataContext::mapWrittenData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasWriteMapping() || not hasMapping());
  if (isMappingRequired()) {
    PRECICE_DEBUG("Map write data \"{}\" from mesh \"{}\"",
                  getDataName(), getMeshName());
    PRECICE_ASSERT(hasMapping());
    _toData->toZero();                                                  // reset _toData
    _mappingContext.mapping->map(_fromData->getID(), _toData->getID()); // map from _fromData to _toData
  }
}

void WriteDataContext::mapWriteDataFrom()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasWriteMapping());
  PRECICE_DEBUG("Map data \"{}\" from mesh \"{}\"", getDataName(), getMeshName());
  _toData->toZero();                                                  // reset _toData
  _mappingContext.mapping->map(_fromData->getID(), _toData->getID()); // store mapped _toData in the _toWaveform
}

} // namespace impl
} // namespace precice
