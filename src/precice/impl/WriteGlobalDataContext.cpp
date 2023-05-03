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

} // namespace precice::impl
