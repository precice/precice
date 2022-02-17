#include "ReadDataContext.hpp"

#include "time/Waveform.hpp"

namespace precice {
namespace impl {

logging::Logger ReadDataContext::_log{"impl::ReadDataContext"};

ReadDataContext::ReadDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    int           interpolationOrder)
    : DataContext(data, mesh)
{
  _providedWaveform = time::PtrWaveform(new time::Waveform(interpolationOrder));
}

void ReadDataContext::configureMapping(const MappingContext &mappingContext, const MeshContext &meshContext)
{
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the read mapping is mapping from needs to be different from _providedData");
  setMapping(mappingContext, data, _providedData);
  PRECICE_ASSERT(hasReadMapping());
}

int ReadDataContext::getInterpolationOrder() const
{
  return _providedWaveform->getInterpolationOrder();
}

void ReadDataContext::storeDataInWaveformFirstSample()
{
  _providedWaveform->storeAtFirstSample(_providedData->values()); // store mapped or received _providedData in the _providedWaveform
}

Eigen::VectorXd ReadDataContext::sampleWaveformAt(double normalizedDt)
{
  return _providedWaveform->sample(normalizedDt);
}

void ReadDataContext::initializeWaveform()
{
  PRECICE_ASSERT(not hasWriteMapping(), "Write mapping does not need waveforms.");
  _providedWaveform->initialize(_providedData->values().size());
  _providedWaveform->storeAtAllSamples(_providedData->values());
}

void ReadDataContext::moveToNextWindow()
{
  _providedWaveform->moveToNextWindow();
}

} // namespace impl
} // namespace precice
