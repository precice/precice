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

int DataContext::getFromDataID() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  PRECICE_ASSERT(_fromData);
  return _fromData->getID();
}

void DataContext::resetProvidedData()
{
  PRECICE_TRACE();
  _providedData->toZero();
  // TODO: consistently reset waveform
  // _providedData->waveform()->toZero();
}

void DataContext::resetToData()
{
  PRECICE_TRACE();
  _toData->toZero();
  // TODO: consistently reset waveform
  // _toData->waveform()->toZero();
}

void DataContext::doWaveformMapping()
{
  resetToData();
  PRECICE_ASSERT(hasMapping());
  PRECICE_DEBUG("Map data \"{}\" from mesh \"{}\"", getDataName(), getMeshName());
  for (int sampleID = 0; sampleID < _providedData->sizeOfSampleStorageInWaveform(); ++sampleID) {
    mapWaveformSample(sampleID);
  }
}

void DataContext::mapWaveformSample(int sampleID)
{
  _fromData->sampleWaveformIntoData(sampleID);                   // put samples from _fromWaveform into _fromData
  mappingContext().mapping->map(getFromDataID(), getToDataID()); // map from _fromData to _toData
  _toData->storeDataInWaveform(sampleID);                        // store _toData at the right place into the _toWaveform
}

int DataContext::getToDataID() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  PRECICE_ASSERT(_toData);
  return _toData->getID();
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

const MappingContext DataContext::mappingContext() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  return _mappingContext;
}

void DataContext::initializeContextWaveforms()
{
  if (hasMapping()) {
    _fromData->initializeWaveform();
    _toData->initializeWaveform();
  } else {
    _providedData->initializeWaveform();
  }
}

void DataContext::sampleWaveformInToData()
{
  PRECICE_TRACE();
  if (hasMapping()) {
    _toData->sampleWaveformIntoData();
  } else {
    _providedData->sampleWaveformIntoData();
  }
}

void DataContext::storeFromDataInWaveform()
{
  PRECICE_TRACE();
  int sampleID = 0;
  if (hasMapping()) {
    _fromData->storeDataInWaveform(sampleID);
  } else {
    _providedData->storeDataInWaveform(sampleID);
  }
}

void DataContext::moveProvidedDataToProvidedWaveformSample(int sampleID)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not hasMapping());
  _providedData->storeDataInWaveform(sampleID);
}

Eigen::VectorXd DataContext::sampleAt(double normalizedDt)
{
  PRECICE_TRACE();
  return _providedData->waveformSampleAt(normalizedDt);
}

void DataContext::moveToNextWindow()
{
  _providedData->moveToNextWindow();
}

} // namespace impl
} // namespace precice
