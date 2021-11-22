#include "precice/impl/DataContext.hpp"
#include <memory>
#include "time/Waveform.hpp"

namespace precice {
namespace impl {

DataContext::DataContext(mesh::PtrData data, mesh::PtrMesh mesh)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(data);
  data->waveform()->initialize(data->values().size());
  data->waveform()->store(data->values());
  _providedData = data;
  PRECICE_ASSERT(data->waveform()->valuesSize() == _providedData->values().size());
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
  PRECICE_ASSERT(_fromData->waveform());
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

int DataContext::getToDataID() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  PRECICE_ASSERT(_toData);
  PRECICE_ASSERT(_toData->waveform());
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
  PRECICE_ASSERT(fromData->waveform() == _providedData->waveform() || toData->waveform() == _providedData->waveform(), "Either _fromData->waveform() or _toData->waveform() has to equal _providedData->waveform().");
  PRECICE_ASSERT(_fromData->waveform()->valuesSize() == _toData->waveform()->valuesSize());
  PRECICE_ASSERT(_toData->waveform() != _fromData->waveform());
}

void DataContext::configureForReadMapping(MappingContext mappingContext, MeshContext fromMeshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(fromMeshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData fromData = fromMeshContext.mesh->data(getDataName());
  PRECICE_ASSERT(fromData != _providedData);
  fromData->waveform()->initialize(fromData->values().size());
  this->setMapping(mappingContext, fromData, _providedData);
  PRECICE_ASSERT(hasReadMapping());
}

void DataContext::configureForWriteMapping(MappingContext mappingContext, MeshContext toMeshContext)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(toMeshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData toData = toMeshContext.mesh->data(getDataName());
  PRECICE_ASSERT(toData != _providedData);
  toData->waveform()->initialize(toData->values().size());
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

void DataContext::initializeProvidedWaveform()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not hasMapping());
  initializeWaveform(_providedData);
}

void DataContext::initializeFromWaveform()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  initializeWaveform(_fromData);
}

void DataContext::initializeToWaveform()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasMapping());
  initializeWaveform(_toData);
}

void DataContext::sampleWaveformInToData()
{
  PRECICE_TRACE();
  if (hasMapping()) {
    sampleWaveformIntoData(_toData);
  } else {
    sampleWaveformIntoData(_providedData);
  }
}

void DataContext::storeFromDataInWaveform()
{
  PRECICE_TRACE();
  if (hasMapping()) {
    storeDataInWaveform(_fromData);
  } else {
    storeDataInWaveform(_providedData);
  }
}

void DataContext::moveWaveformSampleToData(int sampleID)
{
  PRECICE_TRACE();
  sampleWaveformIntoData(_fromData, sampleID);
}

void DataContext::moveDataToWaveformSample(int sampleID)
{
  PRECICE_TRACE();
  storeDataInWaveform(_toData, sampleID);
}

void DataContext::moveProvidedDataToProvidedWaveformSample(int sampleID)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not hasMapping());
  storeDataInWaveform(_providedData, sampleID);
}

int DataContext::sizeOfSampleStorageInWaveform()
{
  PRECICE_TRACE();
  // @todo mpirun -np 4 ./testprecice -t PreciceTests/Serial/MultiCoupling breaks?
  if (hasMapping()) {
    PRECICE_ASSERT(_fromData->waveform()->sizeOfSampleStorage() == _toData->waveform()->sizeOfSampleStorage());
    return _fromData->waveform()->sizeOfSampleStorage();
  } else {
    return _providedData->waveform()->sizeOfSampleStorage();
  }
}

Eigen::VectorXd DataContext::sampleAt(double normalizedDt)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_providedData->waveform()->valuesSize() == _providedData->values().size(),
                 _providedData->waveform()->valuesSize(), _providedData->values().size());

  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");
  return _providedData->waveform()->sample(normalizedDt);
}

void DataContext::initializeWaveform(mesh::PtrData initializingData)
{
  PRECICE_TRACE();
  int sizeOfSampleStorage = sizeOfSampleStorageInWaveform();
  int valuesSize          = initializingData->values().size();
  // PRECICE_ASSERT(valuesSize > 0, valuesSize);  // @todo assertion breaks, but seems like calling advance on empty write data is ok?
  initializingData->waveform()->resizeData(valuesSize);
  for (int sampleID = 0; sampleID < sizeOfSampleStorage; ++sampleID) {
    initializingData->waveform()->storeAt(initializingData->values(), sampleID);
  }
  PRECICE_ASSERT(initializingData->waveform()->valuesSize() == valuesSize);
}

void DataContext::sampleWaveformIntoData(mesh::PtrData data, int sampleID)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(data->waveform()->valuesSize() == data->values().size(),
                 data->waveform()->valuesSize(), data->values().size());
  data->values() = data->waveform()->getSample(sampleID);
}

void DataContext::storeDataInWaveform(mesh::PtrData data, int sampleID)
{
  PRECICE_ASSERT(data->waveform()->valuesSize() == data->values().size(),
                 data->waveform()->valuesSize(), data->values().size());
  data->waveform()->storeAt(data->values(), sampleID);
}

void DataContext::moveProvidedWaveform()
{
  _providedData->waveform()->moveToNextWindow();
}

} // namespace impl
} // namespace precice
