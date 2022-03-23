#include "testing/DataContextFixture.hpp"

namespace precice {
namespace testing {

impl::MappingContext DataContextFixture::mappingContext(precice::impl::DataContext &dataContext)
{
  return dataContext.mappingContext();
}

int DataContextFixture::getProvidedDataID(precice::impl::DataContext &dataContext)
{
  return dataContext.getProvidedDataID();
}

int DataContextFixture::getFromDataID(precice::impl::DataContext &dataContext)
{
  return dataContext.getFromDataID();
}

int DataContextFixture::getToDataID(precice::impl::DataContext &dataContext)
{
    return dataContext.getToDataID();
}

bool DataContextFixture::hasMapping(precice::impl::DataContext &dataContext)
{
  return dataContext.hasMapping();
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
