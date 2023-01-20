#include "ReadDataContext.hpp"

#include "cplscheme/CouplingScheme.hpp"
#include "time/Waveform.hpp"

namespace precice::impl {

logging::Logger ReadDataContext::_log{"impl::ReadDataContext"};

ReadDataContext::ReadDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    int           interpolationOrder)
    : DataContext(data, mesh)
{
  _waveform = std::make_shared<time::Waveform>(interpolationOrder);
}

void ReadDataContext::appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext)
{
  PRECICE_ASSERT(!hasReadMapping(), "The read data context must be unique. Otherwise we would have an ambiguous read data operation on the user side.")
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the read mapping is mapping from needs to be different from _providedData");
  mappingContext.fromData = data;
  mappingContext.toData   = _providedData;
  appendMapping(mappingContext);
  PRECICE_ASSERT(hasReadMapping());
}

int ReadDataContext::getInterpolationOrder() const
{
  return _waveform->getInterpolationOrder();
}

void ReadDataContext::storeDataInWaveform(double relativeDt)
{
  _waveform->store(_providedData->values(), relativeDt); // store mapped or received _providedData in the _waveform
}

Eigen::VectorXd ReadDataContext::sampleWaveformAt(double normalizedDt)
{
  return _waveform->sample(normalizedDt);
}

void ReadDataContext::moveToNextWindow()
{
  _waveform->moveToNextWindow();
}

std::vector<double> ReadDataContext::getReceivedTimes(cplscheme::PtrCouplingScheme cplscheme)
{
  return cplscheme->getReceiveTimes(getDataName());
}

void ReadDataContext::loadReceived(double relativeDt, cplscheme::PtrCouplingScheme cplscheme)
{
  cplscheme->loadReceiveDataFromStorage(getDataName(), relativeDt);
}

} // namespace precice::impl
