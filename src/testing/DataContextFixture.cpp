#include "testing/DataContextFixture.hpp"

namespace precice {
namespace testing {

bool DataContextFixture::hasReadMapping(precice::impl::DataContext &dataContext)
{
  return dataContext.hasReadMapping();
}

bool DataContextFixture::hasWriteMapping(precice::impl::DataContext &dataContext)
{
  return dataContext.hasWriteMapping();
}

void DataContextFixture::resetProvidedData(precice::impl::DataContext &dataContext)
{
  return dataContext.resetProvidedData();
}

void DataContextFixture::resetToData(precice::impl::DataContext &dataContext)
{
  return dataContext.resetToData();
}

} // namespace testing
} // namespace precice
