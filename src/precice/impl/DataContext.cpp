#include "precice/impl/DataContext.hpp"
#include <memory>

namespace precice {
namespace impl {

DataContext::DataContext(mesh::PtrData data, mesh::PtrMesh mesh)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(data);
  _providedData = data;
  PRECICE_ASSERT(mesh);
  _mesh = mesh;
}

mesh::PtrData DataContext::providedData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_providedData);
  return _providedData;
}

std::string DataContext::getDataName() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_providedData);
  return _providedData->getName();
}

int DataContext::getProvidedDataID() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_providedData);
  return _providedData->getID();
}

bool DataContext::isMappingRequired()
{
  using namespace mapping;
  MappingConfiguration::Timing timing    = _mappingContext.timing;
  bool                         hasMapped = _mappingContext.hasMappedData;
  bool                         mapNow    = timing == MappingConfiguration::ON_ADVANCE;
  mapNow |= timing == MappingConfiguration::INITIAL;
  return (hasMapping() && mapNow && (not hasMapped));
}

std::string DataContext::getMeshName() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_mesh);
  return _mesh->getName();
}

int DataContext::getMeshID() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_mesh);
  return _mesh->getID();
}

void DataContext::setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(!hasMapping());
  PRECICE_ASSERT(fromData);
  PRECICE_ASSERT(toData);
  _mappingContext = mappingContext;
  PRECICE_ASSERT(fromData == _providedData || toData == _providedData, "Either fromData or toData has to equal _providedData.");
  PRECICE_ASSERT(fromData->getName() == getDataName());
  _fromData = fromData;
  PRECICE_ASSERT(toData->getName() == getDataName());
  _toData = toData;
  PRECICE_ASSERT(_toData != _fromData);
}

bool DataContext::hasMapping() const
{
  PRECICE_TRACE();
  return hasReadMapping() || hasWriteMapping();
}

bool DataContext::hasReadMapping() const
{
  PRECICE_TRACE();
  return _toData == _providedData;
}

bool DataContext::hasWriteMapping() const
{
  PRECICE_TRACE();
  return _fromData == _providedData;
}

} // namespace impl
} // namespace precice
