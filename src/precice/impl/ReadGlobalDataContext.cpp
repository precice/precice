#include "ReadGlobalDataContext.hpp"

#include "time/Waveform.hpp"

namespace precice::impl {

logging::Logger ReadGlobalDataContext::_log{"impl::ReadGlobalDataContext"};

ReadGlobalDataContext::ReadGlobalDataContext(
    mesh::PtrData data,
    int           interpolationOrder)
    : DataContext(data, nullptr)
{
  _waveform = std::make_shared<time::Waveform>(interpolationOrder);
}

int ReadGlobalDataContext::getInterpolationOrder() const
{
  return _waveform->getInterpolationOrder();
}

void ReadGlobalDataContext::storeDataInWaveform()
{
  _waveform->store(_providedData->values()); // store mapped or received _providedData in the _waveform
}

Eigen::VectorXd ReadGlobalDataContext::sampleWaveformAt(double normalizedDt)
{
  return _waveform->sample(normalizedDt);
}

void ReadGlobalDataContext::initializeWaveform()
{
  PRECICE_ASSERT(not hasWriteMapping(), "Write mapping does not need waveforms.");
  _waveform->initialize(_providedData->values());
}

void ReadGlobalDataContext::moveToNextWindow()
{
  _waveform->moveToNextWindow();
}

} // namespace precice::impl
