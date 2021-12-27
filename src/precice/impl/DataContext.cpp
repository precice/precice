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

void DataContext::mapData(const std::string &mappingType)
{
  using namespace mapping;
  MappingConfiguration::Timing timing    = _mappingContext.timing;
  bool                         hasMapped = _mappingContext.hasMappedData;
  bool                         mapNow    = timing == MappingConfiguration::ON_ADVANCE;
  mapNow |= timing == MappingConfiguration::INITIAL;
  if(hasMapping()){
    if (mapNow && (not hasMapped)) {
      PRECICE_DEBUG("Map \"{}\" data \"{}\" from mesh \"{}\"",
                    mappingType, getDataName(), getMeshName());
      mapWaveformSample();
    }
  } else {
    _providedData->storeDataInWaveform();
  }
}

void DataContext::mapWaveformSample()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  _toData->toZero();                                                  // reset _toData
  _mappingContext.mapping->map(_fromData->getID(), _toData->getID()); // map from _fromData to _toData
  _toData->storeDataInWaveform();                                     // store mapped _toData in the _toWaveform
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
  if (hasMapping()) {
    _fromData->initializeWaveform();
    _toData->initializeWaveform();
  } else {
    _providedData->initializeWaveform();
  }
}

Eigen::VectorXd DataContext::sampleAt(double normalizedDt)
{
  PRECICE_TRACE();
  return _providedData->waveformSampleAt(normalizedDt);
}

void DataContext::moveToNextWindow()
{
  PRECICE_TRACE();
  _providedData->moveToNextWindow();
}

} // namespace impl
} // namespace precice
