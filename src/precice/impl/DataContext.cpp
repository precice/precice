#include "precice/impl/DataContext.hpp"
#include <memory>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "time/Waveform.hpp"

namespace precice {
namespace impl {

DataContext::DataContext(mesh::PtrData data, mesh::PtrMesh mesh)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(data);
  _providedWaveform = data->waveform();
  _providedWaveform->initialize(data->values().size());
  _providedWaveform->store(data->values());
  _providedData = data;
  PRECICE_ASSERT(_providedWaveform->valuesSize() == _providedData->values().size());
  PRECICE_ASSERT(mesh);
  _mesh = mesh;
}

mesh::PtrData DataContext::providedData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_providedData);
  return _providedData;
}

mesh::PtrData DataContext::toData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_toData);
  return _toData;
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
  PRECICE_ASSERT(_fromWaveform);
  return _fromData->getID();
}

void DataContext::resetProvidedData()
{
  PRECICE_TRACE();
  _providedData->toZero();
  // TODO: consistently reset waveform
  // _providedWaveform->toZero();
}

void DataContext::resetToData()
{
  PRECICE_TRACE();
  _toData->toZero();
  // TODO: consistently reset waveform
  // _toWaveform->toZero();
}

int DataContext::getToDataID() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  PRECICE_ASSERT(_toData);
  PRECICE_ASSERT(_toWaveform);
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

void DataContext::setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData, time::PtrWaveform fromWaveform, time::PtrWaveform toWaveform)
{
  PRECICE_TRACE();
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
  PRECICE_ASSERT(_fromWaveform->valuesSize() == _toWaveform->valuesSize());
  PRECICE_ASSERT(_toWaveform != _fromWaveform);
}

void DataContext::configureForReadMapping(MappingContext mappingContext, MeshContext fromMeshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(fromMeshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData fromData = fromMeshContext.mesh->data(getDataName());
  PRECICE_ASSERT(fromData != _providedData);
  time::PtrWaveform ptrFromWaveform = fromData->waveform();
  ptrFromWaveform->initialize(fromData->values().size());
  this->setMapping(mappingContext, fromData, _providedData, ptrFromWaveform, _providedWaveform);
  PRECICE_ASSERT(hasReadMapping());
}

void DataContext::configureForWriteMapping(MappingContext mappingContext, MeshContext toMeshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(toMeshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData toData = toMeshContext.mesh->data(getDataName());
  PRECICE_ASSERT(toData != _providedData);
  time::PtrWaveform ptrToWaveform = toData->waveform();
  ptrToWaveform->initialize(toData->values().size());
  this->setMapping(mappingContext, _providedData, toData, _providedWaveform, ptrToWaveform);
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

void DataContext::initializeProvidedWaveform()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not hasMapping());
  initializeWaveform(_providedData, _providedWaveform);
}

void DataContext::initializeFromWaveform()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  initializeWaveform(_fromData, _fromWaveform);
}

void DataContext::initializeToWaveform()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  initializeWaveform(_toData, _toWaveform);
}

void DataContext::sampleWaveformInToData()
{
  PRECICE_TRACE();
  if (hasMapping()) {
    sampleWaveformIntoData(_toData, _toWaveform);
  } else {
    sampleWaveformIntoData(_providedData, _providedWaveform);
  }
}

void DataContext::storeFromDataInWaveform()
{
  PRECICE_TRACE();
  if (hasMapping()) {
    storeDataInWaveform(_fromData, _fromWaveform);
  } else {
    storeDataInWaveform(_providedData, _providedWaveform);
  }
}

void DataContext::moveWaveformSampleToData(int sampleID)
{
  PRECICE_TRACE();
  sampleWaveformIntoData(_fromData, _fromWaveform, sampleID);
}

void DataContext::moveDataToWaveformSample(int sampleID)
{
  PRECICE_TRACE();
  storeDataInWaveform(_toData, _toWaveform, sampleID);
}

void DataContext::moveProvidedDataToProvidedWaveformSample(int sampleID)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not hasMapping());
  storeDataInWaveform(_providedData, _providedWaveform, sampleID);
}

int DataContext::sizeOfSampleStorageInWaveform()
{
  PRECICE_TRACE();
  // @todo mpirun -np 4 ./testprecice -t PreciceTests/Serial/MultiCoupling breaks?
  if (hasMapping()) {
    PRECICE_ASSERT(_fromWaveform->sizeOfSampleStorage() == _toWaveform->sizeOfSampleStorage());
    return _fromWaveform->sizeOfSampleStorage();
  } else {
    return _providedWaveform->sizeOfSampleStorage();
  }
}

Eigen::VectorXd DataContext::sampleAt(double normalizedDt)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_providedWaveform->valuesSize() == _providedData->values().size(),
                 _providedWaveform->valuesSize(), _providedData->values().size());

  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");
  return _providedWaveform->sample(normalizedDt);
}

void DataContext::initializeWaveform(mesh::PtrData initializingData, time::PtrWaveform initializedWaveform)
{
  PRECICE_TRACE();
  int sizeOfSampleStorage = sizeOfSampleStorageInWaveform();
  int valuesSize          = initializingData->values().size();
  // PRECICE_ASSERT(valuesSize > 0, valuesSize);  // @todo assertion breaks, but seems like calling advance on empty write data is ok?
  initializedWaveform->resizeData(valuesSize);
  for (int sampleID = 0; sampleID < sizeOfSampleStorage; ++sampleID) {
    initializedWaveform->storeAt(initializingData->values(), sampleID);
  }
  PRECICE_ASSERT(initializedWaveform->valuesSize() == valuesSize);
}

void DataContext::sampleWaveformIntoData(mesh::PtrData targetData, time::PtrWaveform sourceWaveform, int sampleID)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(sourceWaveform->valuesSize() == targetData->values().size(),
                 sourceWaveform->valuesSize(), targetData->values().size());
  targetData->values() = sourceWaveform->getSample(sampleID);
}

void DataContext::storeDataInWaveform(mesh::PtrData sourceData, time::PtrWaveform targetWaveform, int sampleID)
{
  PRECICE_ASSERT(targetWaveform->valuesSize() == sourceData->values().size(),
                 targetWaveform->valuesSize(), sourceData->values().size());
  targetWaveform->storeAt(sourceData->values(), sampleID);
}

void DataContext::moveProvidedWaveform()
{
  _providedWaveform->moveToNextWindow();
}

} // namespace impl
} // namespace precice
