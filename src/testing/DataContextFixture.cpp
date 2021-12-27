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

const impl::MappingContext DataContextFixture::mappingContext(precice::impl::DataContext &dataContext) const
{
  return dataContext._mappingContext;
}

int DataContextFixture::getFromDataID(precice::impl::DataContext &dataContext) const
{
  return dataContext._fromData->getID();
}

int DataContextFixture::getToDataID(precice::impl::DataContext &dataContext) const
{
  return dataContext._toData->getID();
}

} // namespace testing
} // namespace precice
