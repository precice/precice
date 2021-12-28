#include "ReadDataContext.hpp"

#include "time/Waveform.hpp"

namespace precice {
namespace impl {
ReadDataContext::ReadDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    int           interpolationOrder)
    : DataContext(data, mesh)
{
  _providedWaveform = time::PtrWaveform(new time::Waveform(interpolationOrder));
}

void ReadDataContext::mapReadData()
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
  _providedWaveform->storeAtFirstSample(_providedData->values()); // store mapped or received _providedData in the _providedWaveform
}

void ReadDataContext::mapReadDataTo()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(hasReadMapping());
  PRECICE_ASSERT(hasMapping());
  PRECICE_DEBUG("Map data \"{}\" to mesh \"{}\"", getDataName(), getMeshName());
  _toData->toZero();                                                  // reset _toData
  _mappingContext.mapping->map(_fromData->getID(), _toData->getID()); // map from _fromData to _toData
  _providedWaveform->storeAtFirstSample(_providedData->values());     // store mapped or received _providedData in the _providedWaveform
}

Eigen::VectorXd ReadDataContext::sampleWaveformAt(double normalizedDt)
{
  PRECICE_TRACE();
  return _providedWaveform->sample(normalizedDt);
}

void ReadDataContext::initializeWaveform()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(not hasWriteMapping(), "Write mapping does not need waveforms.");
  _providedWaveform->initialize(_providedData->values().size());
  _providedWaveform->storeAtAllSamples(_providedData->values());
}

void ReadDataContext::moveToNextWindow()
{
  PRECICE_TRACE();
  _providedWaveform->moveToNextWindow();
}

} // namespace impl
} // namespace precice