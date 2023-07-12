#include "ReadGlobalDataContext.hpp"

namespace precice::impl {

logging::Logger ReadGlobalDataContext::_log{"impl::ReadGlobalDataContext"};

ReadGlobalDataContext::ReadGlobalDataContext(
    mesh::PtrData data)
    : DataContext(data, nullptr)
{
}

void ReadGlobalDataContext::readValue(double normalizedDt, ::precice::span<double> value) const
{
  Eigen::Map<Eigen::MatrixXd>       outputData(value.data(), getDataDimensions(), 1);
  const Eigen::MatrixXd             sample{_providedData->sampleAtTime(normalizedDt)};
  Eigen::Map<const Eigen::MatrixXd> localData(sample.data(), getDataDimensions(), 1);
  outputData.col(0) = localData.col(0);
}

int ReadGlobalDataContext::getWaveformDegree() const
{
  return _providedData->getWaveformDegree();
}

} // namespace precice::impl
