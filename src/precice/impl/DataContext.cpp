#include "precice/impl/DataContext.hpp"
#include <memory>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "time/Waveform.hpp"

namespace precice {
namespace impl {

DataContext::DataContext(mesh::PtrData data, mesh::PtrMesh mesh)
{
  PRECICE_ASSERT(data);
  time::PtrWaveform ptrWaveform(new time::Waveform(data->values().size()));
  _providedWaveform = ptrWaveform;
  _providedWaveform->store(data->values());
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

int DataContext::getFromDataID() const
{
  PRECICE_ASSERT(hasMapping());
  PRECICE_ASSERT(_fromData);
  PRECICE_ASSERT(_fromWaveform);
  return _fromData->getID();
}

void DataContext::resetProvidedData()
{
  _providedData->toZero();
  // TODO: consistently reset waveform
  // _providedWaveform->toZero();
}

void DataContext::resetToData()
{
  _toData->toZero();
  // TODO: consistently reset waveform
  // _toWaveform->toZero();
}

int DataContext::getToDataID() const
{
  PRECICE_ASSERT(hasMapping());
  PRECICE_ASSERT(_toData);
  PRECICE_ASSERT(_toWaveform);
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

void DataContext::setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData, time::PtrWaveform fromWaveform, time::PtrWaveform toWaveform)
{
  PRECICE_ASSERT(!hasMapping());
  PRECICE_ASSERT(fromData);
  PRECICE_ASSERT(toData);
  _mappingContext = mappingContext;
  PRECICE_ASSERT(fromData == _providedData || toData == _providedData, "Either fromData or toData has to equal provided data.");
  PRECICE_ASSERT(fromData->getName() == getDataName());
  _fromData = fromData;
  PRECICE_ASSERT(toData->getName() == getDataName());
  _toData = toData;
  PRECICE_ASSERT(_toData != _fromData);

  PRECICE_ASSERT(fromWaveform);
  PRECICE_ASSERT(toWaveform);
  PRECICE_ASSERT(fromWaveform == _providedWaveform || toWaveform == _providedWaveform, "Either fromWaveform or toWaveform has to equal provided waveform.");
  _fromWaveform = fromWaveform;
  _toWaveform   = toWaveform;
  PRECICE_ASSERT(_toWaveform != _fromWaveform);
}

void DataContext::configureForReadMapping(MappingContext mappingContext, MeshContext meshContext)
{
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData fromData = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(fromData != _providedData);
  time::PtrWaveform ptrFromWaveform(new time::Waveform(fromData->values().size()));
  this->setMapping(mappingContext, fromData, _providedData, ptrFromWaveform, _providedWaveform);
  PRECICE_ASSERT(hasReadMapping());
}

void DataContext::configureForWriteMapping(MappingContext mappingContext, MeshContext meshContext)
{
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData toData = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(toData != _providedData);
  time::PtrWaveform ptrToWaveform(new time::Waveform(toData->values().size()));
  this->setMapping(mappingContext, _providedData, toData, _providedWaveform, ptrToWaveform);
  PRECICE_ASSERT(hasWriteMapping());
}

bool DataContext::hasMapping() const
{
  return hasReadMapping() || hasWriteMapping();
}

bool DataContext::hasReadMapping() const
{
  return _toData == _providedData;
}

bool DataContext::hasWriteMapping() const
{
  return _fromData == _providedData;
}

const MappingContext DataContext::mappingContext() const
{
  PRECICE_ASSERT(hasMapping());
  return _mappingContext;
}

} // namespace impl
} // namespace precice
