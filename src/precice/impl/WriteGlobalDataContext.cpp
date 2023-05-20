#include "WriteGlobalDataContext.hpp"

namespace precice::impl {

logging::Logger WriteGlobalDataContext::_log{"impl::WriteGlobalDataContext"};

WriteGlobalDataContext::WriteGlobalDataContext(
    mesh::PtrData data)
    : DataContext(data, nullptr)
{
}

void WriteGlobalDataContext::resetData(bool atEndOfWindow)
{
  // See also https://github.com/precice/precice/issues/1156.
  _providedData->toZero();

  // reset writeDataBuffer
  _writeDataBuffer.values.setZero();

  if (atEndOfWindow) {
    // need to not only clear _providedData->timeStepsStorage(), but also _toData->timeStepsStorage() as soon as we map from storage to storage.
    _providedData->timeStepsStorage().clear();
  }
}

void WriteGlobalDataContext::writeValueIntoDataBuffer(::precice::span<const double> values)
{
  PRECICE_ASSERT(_writeDataBuffer.values.data());
  Eigen::Map<const Eigen::MatrixXd> inputData(values.data(), getDataDimensions(), 1);
  Eigen::Map<Eigen::MatrixXd>       localData(_writeDataBuffer.values.data(), getDataDimensions(), 1);

  localData.col(0) = inputData.col(0);
}

void WriteGlobalDataContext::storeBufferedData(double currentTime)
{
  _providedData->setSampleAtTime(currentTime, _writeDataBuffer);
}

} // namespace precice::impl
