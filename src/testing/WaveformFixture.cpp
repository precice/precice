#include "testing/WaveformFixture.hpp"

namespace precice {
namespace testing {

int WaveformFixture::numberOfValidSamples(time::Waveform &waveform)
{
  return waveform.numberOfValidSamples();
}

int WaveformFixture::numberOfSamples(time::Waveform &waveform)
{
  return waveform.numberOfSamples();
}

int WaveformFixture::numberOfData(time::Waveform &waveform)
{
  return waveform.numberOfData();
}

double WaveformFixture::getValue(time::Waveform &waveform, int dataID, int sampleID)
{
  return waveform._timeWindows(dataID, sampleID);
}

} // namespace testing
} // namespace precice
