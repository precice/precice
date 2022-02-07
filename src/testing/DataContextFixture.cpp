#include "testing/DataContextFixture.hpp"

namespace precice {
namespace testing {

int DataContextFixture::getProvidedDataID(precice::impl::DataContext &dataContext)
{
  return dataContext._providedData->getID();
}

} // namespace testing
} // namespace precice
