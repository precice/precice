#include "testing/WaveformFixture.hpp"

namespace precice::testing {

int WaveformFixture::numberOfStoredSamples(time::Waveform &waveform)
{
  return waveform._numberOfStoredSamples;
}

int WaveformFixture::maxNumberOfStoredSamples(time::Waveform &waveform)
{
  return waveform.maxNumberOfStoredSamples();
}

int WaveformFixture::valuesSize(time::Waveform &waveform)
{
  return waveform.valuesSize();
}

double WaveformFixture::getValue(time::Waveform &waveform, int valueID, int sampleID)
{
  return waveform._timeWindowsStorage(valueID, sampleID);
}

} // namespace precice::testing
