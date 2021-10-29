#include "testing/WaveformFixture.hpp"

namespace precice {
namespace testing {

int WaveformFixture::numberOfStoredSamples(time::Waveform &waveform)
{
  return waveform._numberOfStoredSamples;
}

int WaveformFixture::sizeOfSampleStorage(time::Waveform &waveform)
{
  return waveform.sizeOfSampleStorage();
}

int WaveformFixture::valuesSize(time::Waveform &waveform)
{
  return waveform.valuesSize();
}

double WaveformFixture::getValue(time::Waveform &waveform, int valueID, int sampleID)
{
  return waveform._timeWindowsStorage(valueID, sampleID);
}

} // namespace testing
} // namespace precice
