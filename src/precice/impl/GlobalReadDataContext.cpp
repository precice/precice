#include "GlobalReadDataContext.hpp"

#include "time/Waveform.hpp"

namespace precice::impl {

logging::Logger GlobalReadDataContext::_log{"impl::GlobalReadDataContext"};

GlobalReadDataContext::GlobalReadDataContext(
    mesh::PtrData data,
    int           interpolationOrder)
    : DataContext(data, nullptr)
{
  _waveform = std::make_shared<time::Waveform>(interpolationOrder);
}

int GlobalReadDataContext::getInterpolationOrder() const
{
  return _waveform->getInterpolationOrder();
}

void GlobalReadDataContext::storeDataInWaveform()
{
  _waveform->store(_providedData->values()); // store mapped or received _providedData in the _waveform
}

Eigen::VectorXd GlobalReadDataContext::sampleWaveformAt(double normalizedDt)
{
  return _waveform->sample(normalizedDt);
}

void GlobalReadDataContext::initializeWaveform()
{
  PRECICE_ASSERT(not hasWriteMapping(), "Write mapping does not need waveforms.");
  _waveform->initialize(_providedData->values());
}

void GlobalReadDataContext::moveToNextWindow()
{
  _waveform->moveToNextWindow();
}

} // namespace precice::impl
