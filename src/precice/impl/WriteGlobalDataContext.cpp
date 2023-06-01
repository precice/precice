#include "WriteGlobalDataContext.hpp"

namespace precice::impl {

logging::Logger WriteGlobalDataContext::_log{"impl::WriteGlobalDataContext"};

WriteGlobalDataContext::WriteGlobalDataContext(
    mesh::PtrData data)
    : DataContext(data, nullptr)
{
  auto dimensions  = getDataDimensions();
  _writeDataBuffer = time::Sample{Eigen::VectorXd(dimensions), Eigen::MatrixXd()};
}

void WriteGlobalDataContext::resetData(bool atEndOfWindow, bool isTimeWindowComplete)
{
  // See also https://github.com/precice/precice/issues/1156.
  _providedData->toZero();

  // reset writeDataBuffer
  _writeDataBuffer.values.setZero();

  if (isTimeWindowComplete) {
    PRECICE_ASSERT(atEndOfWindow, "isTimeWindowComplete without atEndOfWindow is forbidden!");
    auto atEnd = _providedData->timeStepsStorage().stamples().back().sample;
    _providedData->timeStepsStorage().setSampleAtTime(time::Storage::WINDOW_START, atEnd); // manually overwrite value at beginning with value from end. Need this exception for WriteDataContext, because CouplingScheme might not be able to update _providedData, if write mapping sits between _providedData and _toData. CouplingScheme in this case only has access to _toData.
  } else if (atEndOfWindow) {
    _providedData->timeStepsStorage().trim();
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
