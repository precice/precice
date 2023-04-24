#include "GlobalWriteDataContext.hpp"

namespace precice::impl {

logging::Logger GlobalWriteDataContext::_log{"impl::GlobalWriteDataContext"};

GlobalWriteDataContext::GlobalWriteDataContext(
    mesh::PtrData data)
    : DataContext(data, nullptr)
{
}

mesh::PtrData GlobalWriteDataContext::providedData()
{
  PRECICE_ASSERT(_providedData);
  return _providedData;
}

} // namespace precice::impl
