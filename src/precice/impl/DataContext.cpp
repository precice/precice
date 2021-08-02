#include "precice/impl/DataContext.hpp"
#include <memory>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace impl {

DataContext::DataContext(mesh::PtrData data, mesh::PtrMesh mesh)
{
  PRECICE_ASSERT(data);
  _providedData = data;
  PRECICE_ASSERT(mesh);
  _mesh = mesh;
}

mesh::PtrData DataContext::providedData()
{
  PRECICE_ASSERT(_providedData);
  return _providedData;
}

std::string DataContext::getDataName() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->getName();
}

int DataContext::getProvidedDataID() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->getID();
}

mesh::PtrData DataContext::fromData()
{
  PRECICE_ASSERT(_hasMapping);
  PRECICE_ASSERT(_fromData);
  return _fromData;
}

int DataContext::getFromDataID() const
{
  PRECICE_ASSERT(_hasMapping);
  PRECICE_ASSERT(_fromData);
  return _fromData->getID();
}

mesh::PtrData DataContext::toData()
{
  PRECICE_ASSERT(_hasMapping);
  PRECICE_ASSERT(_toData);
  return _toData;
}

int DataContext::getToDataID() const
{
  PRECICE_ASSERT(_hasMapping);
  PRECICE_ASSERT(_toData);
  return _toData->getID();
}

std::string DataContext::getMeshName() const
{
  PRECICE_ASSERT(_mesh);
  return _mesh->getName();
}

int DataContext::getMeshID() const
{
  PRECICE_ASSERT(_mesh);
  return _mesh->getID();
}

void DataContext::setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData)
{
  PRECICE_ASSERT(!_hasMapping);
  _hasMapping = true;
  PRECICE_ASSERT(fromData);
  PRECICE_ASSERT(toData);
  _mappingContext = mappingContext;
  PRECICE_ASSERT(fromData == _providedData || toData == _providedData, "Either fromData or toData has to equal provided data.");
  PRECICE_ASSERT(fromData->getName() == getDataName());
  _fromData = fromData;
  PRECICE_ASSERT(toData->getName() == getDataName());
  _toData = toData;
}

void DataContext::configureForReadMapping(MappingContext mappingContext, MeshContext meshContext)
{
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData fromData = meshContext.mesh->data(getDataName());
  this->setMapping(mappingContext, fromData, _providedData);
}

void DataContext::configureForWriteMapping(MappingContext mappingContext, MeshContext meshContext)
{
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData toData = meshContext.mesh->data(getDataName());
  this->setMapping(mappingContext, _providedData, toData);
}

bool DataContext::hasMapping() const
{
  return _hasMapping;
}

const MappingContext DataContext::mappingContext() const
{
  PRECICE_ASSERT(_hasMapping);
  return _mappingContext;
}

} // namespace impl
} // namespace precice
