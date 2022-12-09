#include "testing/DataContextFixture.hpp"

namespace precice::testing {

std::vector<impl::MappingContext> DataContextFixture::mappingContexts(precice::impl::DataContext &dataContext)
{
  return dataContext._mappingContexts;
}

int DataContextFixture::getProvidedDataID(precice::impl::DataContext &dataContext)
{
  return dataContext._providedData->getID();
}

int DataContextFixture::getFromDataID(precice::impl::DataContext &dataContext, int dataVectorIndex)
{
  return dataContext._mappingContexts[dataVectorIndex].fromData->getID();
}

int DataContextFixture::getToDataID(precice::impl::DataContext &dataContext, int dataVectorIndex)
{
  return dataContext._mappingContexts[dataVectorIndex].toData->getID();
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

} // namespace precice::testing
