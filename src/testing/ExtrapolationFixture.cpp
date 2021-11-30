#include "testing/ExtrapolationFixture.hpp"

namespace precice {
namespace testing {

int ExtrapolationFixture::numberOfStoredSamples(cplscheme::Extrapolation &extrapolation)
{
  return extrapolation._numberOfStoredSamples;
}

int ExtrapolationFixture::sizeOfSampleStorage(cplscheme::Extrapolation &extrapolation)
{
  return extrapolation.sizeOfSampleStorage();
}

int ExtrapolationFixture::valuesSize(cplscheme::Extrapolation &extrapolation)
{
  return extrapolation.valuesSize();
}

double ExtrapolationFixture::getValue(cplscheme::Extrapolation &extrapolation, int valueID, int sampleID)
{
  return extrapolation._timeWindowsStorage(valueID, sampleID);
}

} // namespace testing
} // namespace precice
