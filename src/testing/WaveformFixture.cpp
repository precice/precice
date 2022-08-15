#include "testing/WaveformFixture.hpp"

namespace precice {
namespace testing {

int WaveformFixture::numberOfStoredSamples(time::Waveform &waveform)
{
  return waveform._numberOfStoredSamples;
}

int WaveformFixture::maxNumberOfStoredWindows(time::Waveform &waveform)
{
  return waveform.maxNumberOfStoredWindows();
}

int WaveformFixture::valuesSize(time::Waveform &waveform)
{
  return waveform.valuesSize();
}

} // namespace testing
} // namespace precice
