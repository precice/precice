#include "ReadGlobalDataContext.hpp"

#include "time/Waveform.hpp"

namespace precice::impl {

logging::Logger ReadGlobalDataContext::_log{"impl::ReadGlobalDataContext"};

ReadGlobalDataContext::ReadGlobalDataContext(
    mesh::PtrData data,
    int           interpolationOrder)
    : DataContext(data, nullptr)
{
  _waveform = std::make_shared<time::Waveform>(interpolationOrder, data);
}

void ReadGlobalDataContext::readValue(double normalizedDt, ::precice::span<double> value) const
{
  Eigen::Map<Eigen::MatrixXd>       outputData(value.data(), getDataDimensions(), 1);
  const Eigen::MatrixXd             sample{_waveform->sample(normalizedDt)};
  Eigen::Map<const Eigen::MatrixXd> localData(sample.data(), getDataDimensions(), 1);
  outputData.col(0) = localData.col(0);
}

int ReadGlobalDataContext::getInterpolationOrder() const
{
  return _waveform->getInterpolationOrder();
}

Eigen::VectorXd ReadGlobalDataContext::sampleWaveformAt(double normalizedDt)
{
  return _waveform->sample(normalizedDt);
}

void ReadGlobalDataContext::moveToNextWindow()
{
  _waveform->moveToNextWindow();
}

} // namespace precice::impl
