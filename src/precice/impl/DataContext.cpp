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

void DataContext::mapWrittenData()
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

void DataContext::mapReadData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasReadMapping() || not hasMapping());
  if (isMappingRequired()) {
    PRECICE_DEBUG("Map read data \"{}\" from mesh \"{}\"",
                  getDataName(), getMeshName());
    PRECICE_ASSERT(hasMapping());
    _toData->toZero();                                                  // reset _toData
    _mappingContext.mapping->map(_fromData->getID(), _toData->getID()); // map from _fromData to _toData
  }
  _providedData->storeDataInWaveform(); // store mapped or received _providedData in the _providedWaveform
}

void DataContext::mapWriteDataFrom()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasWriteMapping());
  PRECICE_DEBUG("Map data \"{}\" from mesh \"{}\"", getDataName(), getMeshName());
  _toData->toZero();                                                  // reset _toData
  _mappingContext.mapping->map(_fromData->getID(), _toData->getID()); // store mapped _toData in the _toWaveform
}

void DataContext::mapReadDataTo()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasReadMapping());
  PRECICE_ASSERT(hasMapping());
  PRECICE_DEBUG("Map data \"{}\" to mesh \"{}\"", getDataName(), getMeshName());
  _toData->toZero();                                                  // reset _toData
  _mappingContext.mapping->map(_fromData->getID(), _toData->getID()); // map from _fromData to _toData
  _providedData->storeDataInWaveform();                               // store mapped _toData in the _toWaveform
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

void DataContext::configureForReadMapping(MappingContext mappingContext, MeshContext fromMeshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(fromMeshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData fromData = fromMeshContext.mesh->data(getDataName());
  PRECICE_ASSERT(fromData != _providedData);
  this->setMapping(mappingContext, fromData, _providedData);
  PRECICE_ASSERT(hasReadMapping());
}

void DataContext::configureForWriteMapping(MappingContext mappingContext, MeshContext toMeshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(toMeshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData toData = toMeshContext.mesh->data(getDataName());
  PRECICE_ASSERT(toData != _providedData);
  this->setMapping(mappingContext, _providedData, toData);
  PRECICE_ASSERT(hasWriteMapping());
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

void DataContext::initializeContextWaveforms()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not hasWriteMapping(), "Write mapping does not need waveforms.");
  _providedData->initializeWaveform();
}

void DataContext::moveToNextWindow()
{
  PRECICE_TRACE();
  _providedData->moveToNextWindow();
}

} // namespace impl
} // namespace precice
