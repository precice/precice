#include "precice/impl/DataContext.hpp"
#include <memory>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace impl {

DataContext::DataContext(mesh::PtrData data, mesh::PtrMesh mesh)
{
  PRECICE_ASSERT(data);
  _participantData = data;
  PRECICE_ASSERT(mesh);
  _mesh = mesh;
}

mesh::PtrData DataContext::participantData()
{
  PRECICE_ASSERT(_participantData);
  return _participantData;
}

std::string DataContext::getDataName() const
{
  PRECICE_ASSERT(_participantData);
  return _participantData->getName();
}

int DataContext::getParticipantDataID() const
{
  PRECICE_ASSERT(_participantData);
  return _participantData->getID();
}

mesh::PtrData DataContext::fromData()
{
  PRECICE_ASSERT(_fromData);
  return _fromData;
}

int DataContext::getFromDataID() const
{
  PRECICE_ASSERT(_fromData);
  return _fromData->getID();
}

mesh::PtrData DataContext::toData()
{
  PRECICE_ASSERT(_toData);
  return _toData;
}

int DataContext::getToDataID() const
{
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
  PRECICE_ASSERT(!hasMapping());
  PRECICE_ASSERT(fromData);
  PRECICE_ASSERT(toData);
  _mappingContext = mappingContext;
  PRECICE_ASSERT(fromData->getName() == getDataName());
  _fromData       = fromData;
  PRECICE_ASSERT(toData->getName() == getDataName());
  _toData         = toData;
}

bool DataContext::hasMapping() const
{
  return _toData && _fromData && (_toData != _fromData);
}

const MappingContext DataContext::mappingContext() const
{
  return _mappingContext;
}

} // namespace impl
} // namespace precice
