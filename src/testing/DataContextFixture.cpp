#include "testing/DataContextFixture.hpp"

namespace precice {
namespace testing {

std::vector<impl::MappingContext> DataContextFixture::mappingContext(precice::impl::DataContext &dataContext)
{
  return dataContext._mappingContexts;
}

int DataContextFixture::getProvidedDataID(precice::impl::DataContext &dataContext)
{
  return dataContext._providedData->getID();
}

int DataContextFixture::getFromDataID(precice::impl::DataContext &dataContext, int dataVectorIndex)
{
  return dataContext.getFromDataID(dataVectorIndex);
}

int DataContextFixture::getToDataID(precice::impl::DataContext &dataContext, int dataVectorIndex)
{
  return dataContext.getToDataID(dataVectorIndex);
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
