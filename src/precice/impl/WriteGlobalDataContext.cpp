#include "WriteGlobalDataContext.hpp"

namespace precice::impl {

logging::Logger WriteGlobalDataContext::_log{"impl::WriteGlobalDataContext"};

WriteGlobalDataContext::WriteGlobalDataContext(
    mesh::PtrData data)
    : DataContext(data, nullptr)
{
}

mesh::PtrData WriteGlobalDataContext::providedData()
{
  PRECICE_ASSERT(_providedData);
  return _providedData;
}

void WriteGlobalDataContext::writeValue(::precice::span<const double> values)
{
  Eigen::Map<const Eigen::MatrixXd> inputData(values.data(), getDataDimensions(), 1);
  Eigen::Map<Eigen::MatrixXd>       localData(_providedData->values().data(), getDataDimensions(), 1);

  localData.col(0) = inputData.col(0);
}

} // namespace precice::impl
