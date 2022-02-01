#include "testing/DataContextFixture.hpp"

namespace precice {
namespace testing {

int DataContextFixture::getProvidedDataID(precice::impl::DataContext &dataContext)
{
  return dataContext.getProvidedDataID();
}

bool DataContextFixture::hasReadMapping(precice::impl::DataContext &dataContext)
{
  return dataContext.hasReadMapping();
}

bool DataContextFixture::hasWriteMapping(precice::impl::DataContext &dataContext)
{
  return dataContext.hasWriteMapping();
}

} // namespace testing
} // namespace precice
