#include "testing/ExtrapolationFixture.hpp"

namespace precice {
namespace testing {

int ExtrapolationFixture::numberOfStoredSamples(cplscheme::impl::Extrapolation &extrapolation)
{
  return extrapolation._numberOfStoredSamples;
}

int ExtrapolationFixture::sizeOfSampleStorage(cplscheme::impl::Extrapolation &extrapolation)
{
  return extrapolation.sizeOfSampleStorage();
}

int ExtrapolationFixture::valuesSize(cplscheme::impl::Extrapolation &extrapolation)
{
  return extrapolation.valuesSize();
}

double ExtrapolationFixture::getValue(cplscheme::impl::Extrapolation &extrapolation, int valueID, int sampleID)
{
  return extrapolation._timeWindowsStorage(valueID, sampleID);
}

} // namespace testing
} // namespace precice
